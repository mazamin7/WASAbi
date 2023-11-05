clear all, close all, clc;

addpath data\
addpath functions\


% Parameters
source_filename = 'source_0.txt';
receiver_filename = 'response_0.txt';
room_filename = 'hall.txt';

c = 343.5;

sim_dur = 2;
% sim_dur = 5;

dh = 0.2;
% dh = 0.05;
% dh = 0.1;

% boundary absorption coefficient
alpha_b = 1.0;
% alpha_b = 0.9;
% alpha_b = 0.8;
% alpha_b = 0.5;
% alpha_b = 0.2;

% air absorption coefficient
% alpha_a = 10;
alpha_a = 0;

fs_new = 44100;

% End parameters


switch dh
    case 0.05
        fs = 16000; % Sampling rate of impulse response
        fcut = 3000; % Bandwidth of impulse response
    case 0.1
        fs = 8000;
        fcut = 1500;
    case 0.2
        fs = 5000;
        fcut = 750;
    case 0.5
        fs = 1600;
        fcut = 300;
end

% Compute S and V
[S, V] = compute_geometry_S_V(room_filename)

% Load source and response
src = load(source_filename);
rr = load(receiver_filename);
Ns = size(rr,1);

% Truncate response before there is only noise
A = S * alpha_b;

% Sabine's formula
T60_sabine = 0.161 * V / (A + 4*alpha_a*V);
T60_sabine_sample = floor(T60_sabine * fs);
disp(['Sabine''s T60: ' num2str(T60_sabine) ' sec.']);
% rr(T60_sabine_sample+1:end) = 0;

% Eyring's formula
T60_eyring = 0.161 * V / (S * abs(log(1 - alpha_b)) + 4*alpha_a*V);
T60_eyring_sample = floor(T60_eyring * fs);
disp(['Eyring''s T60: ' num2str(T60_eyring) ' sec.']);
% rr(T60_eyring_sample+1:end) = 0;

% Compute impulse response through deconvolution
rir = xcorr(rr,src(end:-1:1));
rir = rir(1:(end+1)/2);
b = fir1(1024,fcut/(fs/2));
rir = fftfilt(b,rir);
rir = normalizeIR(rir);

% Upsample
imp_up = resample(rir, fs_new, fs);  % Upsample imp1 to fs2
imp_up = normalizeIR(imp_up);

% Perform peak detection
[~, peakIndices] = findpeaks(imp_up, fs_new, 'MinPeakWidth', 0.001);

% Extract amplitude and delay information from the peaks
peakAmplitudes = imp_up(floor(peakIndices*fs_new));
peakDelays = (peakIndices - 1) / fs_new;

% Create upsampled impulse response containing only the peaks
imp_up_clean = zeros(size(imp_up));
imp_up_clean(floor(peakIndices*fs_new)) = imp_up(floor(peakIndices*fs_new));

% High-pass filtering
[b, a] = butter(2, fcut / (fs_new/2), 'high');  % High-pass filter coefficients
imp_up_clean_high = filter(b, a, imp_up_clean);
imp_up_clean_high = normalizeIR(imp_up_clean_high);

% Sum low frequency impulse response and high frequency impulse response
imp_sum = imp_up + imp_up_clean_high * 10;
rir2 = normalizeIR(imp_sum);

t_axis = (0:length(rir)-1) / fs;


% Plotting the result
fig = figure();

subplot(211);
plot(t_axis, rir);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 sim_dur]);
ylim([-1,1]);
% ylim([-80,0]);
title('Room''s point-to-point impulse response');

% Compute FFT
N = length(rir);  % Length of the signals for FFT
f_axis = (-N/2:N/2-1) * (fs/N);  % Frequency axis

% Compute FFT of rir
rir_fft = abs(fftshift(fft(rir, N)));

% Plot FFT
figure(fig);
subplot(212);
loglog(f_axis, rir_fft);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0, fcut]);
% xticks([0 1 2 3 4 5] * 500);
ylim([1e-2, 1e3]);
yticks([1e-2 1e-1 1 1e1 1e2 1e3]);
title('Room''s point-to-point frequency response');


% Energy decay curve

% Calculate the squared impulse response
h_squared = rir.^2;

% Compute the cumulative integral from 0 to t1(end)
cum_integral = cumtrapz(t_axis, h_squared);

% Compute the EDC by subtracting the cumulative integral at t1(end) from the total cumulative integral
EDC = cum_integral(end) - cum_integral;


% Theoretical decay
max_peak = EDC(1);
% t0_sample = find(EDC < max_peak*0.99, 1);
% t0 = t0_sample / fs;
t0_sample = find(rr == 0, 1);
t0 = t0_sample/fs;
decay_theory = max_peak * exp(-c*A/4/V*(t_axis - t0)) .* exp(-2*(t_axis-t0)*alpha_a) .* (t_axis >= t0) + max_peak * (t_axis < t0);


% Plotting EDC
fig_EDC = figure();

subplot(211);
plot(t_axis, EDC, 'DisplayName', 'Numerical');
hold on;
plot(t_axis, decay_theory, 'r--', 'DisplayName', 'Theory');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 sim_dur]);
% ylim([-1,1]);
% ylim([-80,0]);
title('Room''s point-to-point energy decay curve');
legend('show');

subplot(212);
plot(t_axis, 10*log10(EDC) - 10*log10(EDC(1)), 'DisplayName', 'Numerical');
hold on;
plot(t_axis, 10*log10(decay_theory) - 10*log10(EDC(1)), 'r--', 'DisplayName', 'Theory');
xlabel('Time (s)');
ylabel('Amplitude (dB)');
xlim([0 sim_dur]);
% ylim([-1,1]);
ylim([-80,0]);
title('Room''s point-to-point energy decay curve in dB');
legend('show');


% Plotting the upsampled result
t_axis2 = (0:length(rir2)-1) / fs_new;

fig2 = figure();

subplot(211);
plot(t_axis2, rir2);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 sim_dur]);
ylim([-1,1]);
% ylim([-80,0]);
title('Room''s point-to-point upsampled impulse response');

% Compute FFT
N = length(rir2);  % Length of the signals for FFT
f_axis2 = (-N/2:N/2-1) * (fs_new/N);  % Frequency axis

% Compute FFT of rir
rir2_fft = abs(fftshift(fft(rir2, N)));

% Plot FFT
figure(fig2);
subplot(212);
loglog(f_axis2, rir2_fft);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0, fs_new/2]);
% xticks([0 1 2 3 4 5] * 500);
ylim([1e-2, 1e3]);
yticks([1e-2 1e-1 1 1e1 1e2 1e3]);
title('Room''s point-to-point upsampled frequency response');


% Computing T30 

% t1: source turned off
t1_sample = find(rr == 0, 1);
t1 = t1_sample/fs;

% t2: db(EDC) = -30 dB
t2_sample = find(abs(10*log10(EDC) - 10*log10(EDC(1)) - (-30)) < 0.5, 1);
t2 = t2_sample/fs;

T30 = 2 * (t2 - t1);
disp(['Computed T30: ' num2str(T30) ' sec.']);

% figure(fig_EDC);
% subplot(212);
% hold on;
% % Add a vertical line at t30
% line([T30+t1 T30+t1], ylim, 'Color', 'r', 'LineStyle', '--');


% Computing T60 

% t1: source turned off
t1_sample = find(rr == 0, 1);
t1 = t1_sample/fs;

% t2: db(EDC) = -60 dB
t2_sample = find(abs(10*log10(EDC) - 10*log10(EDC(1)) - (-60)) < 0.5, 1);
t2 = t2_sample/fs;

T60 = t2 - t1;
disp(['Computed T60: ' num2str(T60) ' sec.']);

% figure(fig_EDC);
% subplot(212);
% hold on;
% % Add a vertical line at t30
% line([T60+t1 T60+t1], ylim, 'Color', 'r', 'LineStyle', '--');


% Auralization
[x, ori_sr] = audioread('welcome.wav');

x = resample(x, fs_new, ori_sr);
y = fftfilt(rir2, x);
y = normalizeIR(y);

audiowrite('welcome_reverb.wav', y, fs_new);






