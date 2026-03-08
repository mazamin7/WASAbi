clear all, close all, clc;

% Default experiment path
if ~exist('experiment_path', 'var')
    experiment_path = '../../experiments/hall';
end

source_filename = [experiment_path, '/output/source_data_0.bin'];

% Load Experiment Config
config_str = fileread([experiment_path, '/config.json']);
config = jsondecode(config_str);

% Parameters
c = config.simulation.c0;
sim_dur = config.simulation.duration;
dh = config.simulation.dh;
dt = config.simulation.dt;

fcut = 150 / dh;
fs = round(1/dt); % Sampling rate of impulse response

% Load binary source
fid_src = fopen(source_filename, 'r');
src = fread(fid_src, inf, 'double');
fclose(fid_src);
Ns = size(src,1);

t_axis = (0:length(src)-1) / fs;


% Plotting the result
fig = figure();

subplot(211);
plot(t_axis, src);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 sim_dur]);
% ylim([-1,1]);
% ylim([-80,0]);
title('Impulse');


% Compute FFT
N = length(src);  % Length of the signals for FFT
f_axis = (-N/2:N/2-1) * (fs/N);  % Frequency axis

% Compute FFT of src
rir_fft = abs(fftshift(fft(src, N)));

% Plot FFT
figure(fig);
subplot(212);
loglog(f_axis, rir_fft);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0, fs/2]);
% xticks([0 1 2 3 4 5] * 500);
% ylim([1e-2, 1e3]);
% yticks([1e-2 1e-1 1 1e1 1e2 1e3]);
title('FFT of the impulse');

% Add vertical dashed red line at fcut
hold on;
plot([fcut, fcut], ylim, 'r--');
hold off;




