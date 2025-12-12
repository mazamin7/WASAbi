clear all; close all; clc;

source_filename = 'source_0.txt';
receiver_filename = 'response_0.txt';
alpha1 = 0; % air telegrapher damping
alpha2 = 0; % air viscoelastic damping
alpha_b = 1; % boundary absorption

% Parameters
c = 343.5;
sim_dur = 1;
dh = 0.5;


% Sampling rate and bandwidth based on dh
switch dh
    case 0.05
        fs = 16000; % Sampling rate of impulse response
    case 0.1
        fs = 8000;
    case 0.2
        fs = 5000;
    case 0.5
        fs = 1600;
end

% Compute S and V
S = 2720; V = 7600;

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

    % Load source and response
    src = load(source_filename);
    rr = load(receiver_filename);
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