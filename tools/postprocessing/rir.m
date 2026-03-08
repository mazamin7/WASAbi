% clear all;
close all; clc;

% Default experiment path
if ~exist('experiment_path', 'var')
    experiment_path = '../../experiments/hall';
end

% Correct Binary Filenames
source_filename = [experiment_path, '/output/source_data_0.bin'];
receiver_filename = [experiment_path, '/output/response_data_0.bin'];

% Load Asset Config
asset_str = fileread([experiment_path, '/asset.json']);
asset = jsondecode(asset_str);

% Load Experiment Config
config_str = fileread([experiment_path, '/config.json']);
config = jsondecode(config_str);

alpha1 = asset.medium_properties.alpha1;
alpha2 = asset.medium_properties.alpha2;
c = asset.medium_properties.c0;
sim_dur = config.simulation.duration;
dh = config.simulation.dh;
dt = config.simulation.dt;

% Calculate Surface Area-Weighted alpha_b from all partitions
total_area = 0;
weighted_alpha = 0;
V = 0;

for p_idx = 1:length(asset.partitions)
    p = asset.partitions(p_idx);
    
    % Raw surface areas of faces (w=x, h=y, d=z dimensions in meters)
    area_x = p.h * p.d;
    area_y = p.w * p.d;
    area_z = p.w * p.h;
    
    % Calculate overlaps with all other partitions to exclude air-to-air interfaces
    overlap_x_minus = 0; overlap_x_plus = 0;
    overlap_y_minus = 0; overlap_y_plus = 0;
    overlap_z_minus = 0; overlap_z_plus = 0;
    
    for o_idx = 1:length(asset.partitions)
        if o_idx == p_idx
            continue;
        end
        o = asset.partitions(o_idx);
        
        % Overlap on X faces
        if p.x == o.x + o.w
            intersect_y = max(0, min(p.y + p.h, o.y + o.h) - max(p.y, o.y));
            intersect_z = max(0, min(p.z + p.d, o.z + o.d) - max(p.z, o.z));
            overlap_x_minus = overlap_x_minus + (intersect_y * intersect_z);
        end
        if p.x + p.w == o.x
            intersect_y = max(0, min(p.y + p.h, o.y + o.h) - max(p.y, o.y));
            intersect_z = max(0, min(p.z + p.d, o.z + o.d) - max(p.z, o.z));
            overlap_x_plus = overlap_x_plus + (intersect_y * intersect_z);
        end
        
        % Overlap on Y faces
        if p.y == o.y + o.h
            intersect_x = max(0, min(p.x + p.w, o.x + o.w) - max(p.x, o.x));
            intersect_z = max(0, min(p.z + p.d, o.z + o.d) - max(p.z, o.z));
            overlap_y_minus = overlap_y_minus + (intersect_x * intersect_z);
        end
        if p.y + p.h == o.y
            intersect_x = max(0, min(p.x + p.w, o.x + o.w) - max(p.x, o.x));
            intersect_z = max(0, min(p.z + p.d, o.z + o.d) - max(p.z, o.z));
            overlap_y_plus = overlap_y_plus + (intersect_x * intersect_z);
        end
        
        % Overlap on Z faces
        if p.z == o.z + o.d
            intersect_x = max(0, min(p.x + p.w, o.x + o.w) - max(p.x, o.x));
            intersect_y = max(0, min(p.y + p.h, o.y + o.h) - max(p.y, o.y));
            overlap_z_minus = overlap_z_minus + (intersect_x * intersect_y);
        end
        if p.z + p.d == o.z
            intersect_x = max(0, min(p.x + p.w, o.x + o.w) - max(p.x, o.x));
            intersect_y = max(0, min(p.y + p.h, o.y + o.h) - max(p.y, o.y));
            overlap_z_plus = overlap_z_plus + (intersect_x * intersect_y);
        end
    end
    
    % Net areas for absorption weighting (interfaces have 0 absorption contribution here)
    net_x_minus = area_x - overlap_x_minus;
    net_x_plus = area_x - overlap_x_plus;
    net_y_minus = area_y - overlap_y_minus;
    net_y_plus = area_y - overlap_y_plus;
    net_z_minus = area_z - overlap_z_minus;
    net_z_plus = area_z - overlap_z_plus;

    % Volume
    V = V + (p.w * p.h * p.d);
    
    % Aggregate absorption (MUST be explicit)
    ba = p.boundary_absorption;
    weighted_alpha = weighted_alpha + ba.x_minus * net_x_minus;
    weighted_alpha = weighted_alpha + ba.x_plus * net_x_plus;
    weighted_alpha = weighted_alpha + ba.y_minus * net_y_minus;
    weighted_alpha = weighted_alpha + ba.y_plus * net_y_plus;
    weighted_alpha = weighted_alpha + ba.z_minus * net_z_minus;
    weighted_alpha = weighted_alpha + ba.z_plus * net_z_plus;
    
    total_area = total_area + net_x_minus + net_x_plus + net_y_minus + net_y_plus + net_z_minus + net_z_plus;
end

alpha_b = weighted_alpha / total_area;
S = total_area;

% Sampling rate based on dt
fs = round(1 / dt);

% Compute equivalent absorption area
A = S * alpha_b;

% Define octave bands (central frequencies)
octave_bands = [125, 250, 500, 1000, 2000, 4000];

for i = 1:length(octave_bands)
    fc = octave_bands(i);
    f_low = fc / sqrt(2);
    f_high = fc * sqrt(2);
    
    % Check if octave band is within valid frequency range
    if f_high >= fs/2
        disp(['Skipping octave band ' num2str(fc) ' Hz because it exceeds Nyquist frequency.']);
        continue;
    end

    % Create a figure for the current octave band
    figure_edc = figure;
    
    alpha_a = alpha1 + alpha2 * (2*pi*fc)^2;
    disp(['Air absorption: ' num2str(alpha_a) ' 1/sec.']);

    % Load Binary Source
    fid_src = fopen(source_filename, 'r');
    src = fread(fid_src, inf, 'double');
    fclose(fid_src);

    % Stream Response Binary
    fileID = fopen(receiver_filename, 'r');
    rr = fread(fileID, inf, 'double');
    fclose(fileID);

    Ns = size(rr,1);
    
    % Compute impulse response through deconvolution
    rir_ = xcorr(rr,src(end:-1:1));
    rir_ = rir_(1:(end+1)/2);
    rir_ = normalizeIR(rir_);
    
    t_axis = (0:length(rir_)-1) / fs;
    
    % Compute FFT of unfiltered impulse response
    N = length(rir_);  % Length of the signals for FFT
    f_axis = (-N/2:N/2-1) * (fs/N);  % Frequency axis
    rir_fft = abs(fftshift(fft(rir_, N)));
        
    % Filter the impulse response
    [b, a] = butter(4, [f_low, f_high] / (fs/2));
    rir_filt = filtfilt(b, a, rir_);

    % Energy decay curve
    h_squared = rir_filt.^2;
    cum_integral = cumtrapz(t_axis, h_squared);
    EDC = cum_integral(end) - cum_integral;

    % Theoretical decay
    max_peak = EDC(1);
    t0_sample = find(EDC < max_peak*0.99, 1);
    t0 = t0_sample / fs;
    decay_theory = max_peak * exp(-c*A/4/V*(t_axis - t0)) .* exp(-2*(t_axis-t0)*alpha_a) .* (t_axis >= t0) + max_peak * (t_axis < t0);

    % Plotting EDC
    figure(figure_edc);
    plot(t_axis, EDC, 'DisplayName', 'Numerical', 'LineWidth', 2);
    hold on;
    plot(t_axis, decay_theory, 'r--', 'DisplayName', 'Theory', 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Amplitude', 'FontSize', 14);
    xlim([0 sim_dur]);
    title(['Octave Band ', num2str(fc), ' Hz'], 'FontSize', 14);
    legend('show', 'FontSize', 14);
    
    figure(figure_edc);
    figure_edc.Position = [100 100 600 400];
    
    disp("------------------------");
end