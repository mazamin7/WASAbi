clear all, close all, clc;

addpath data\
addpath functions\

source_filename = 'source_0.txt';

% Parameters
c = 343.5;

sim_dur = 8e-3;
% sim_dur = 5;

dh = 0.2;
% dh = 0.05;
% dh = 0.1;

switch dh
    case 0.05
        dt = 0.625e-4;
        fcut = 3000; % Bandwidth of impulse response
    case 0.1
        dt = 1.25e-4;
        fcut = 1500;
    case 0.2
        dt = 2e-4;
        fcut = 750;
    case 0.5
        dt = 6.25e-4;
        fcut = 300;
end

fs = 1/dt; % Sampling rate of impulse response

% Load source
src = load(source_filename);
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




