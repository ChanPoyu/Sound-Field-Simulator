close all;
clear all;

%% Gaussian band noise

N = 5000;

fs = 48000;

f = linspace(0, 24000, round(10327/ 2));

Signal = zeros(1, N * 2);

center_frequency = 4000;

center_sample = 2500 * (center_frequency / fs);

amp_max = 1;

a = 50;

for n = 1 : 1 : N
   r_sqr = (n - center_sample) ^ 2;
   Signal_amp = amp_max * exp(-r_sqr / a^3);
%    Signal(n) = Signal_amp * exp (1i * 4 * pi * ((n- center_sample) ^ 2) / (N ^ 2));
   Signal(n) = Signal_amp * exp (-1i * 4 * pi * n^2 / N);
   
end



S_temp = flip(Signal(1:N));
S_temp = conj(S_temp);

Signal(N+1 : end) = S_temp;


%% LPF
LPF = dsp.LowpassFilter;
LPF.SampleRate = 48000;
LPF.StopbandFrequency = 4000;
LPF.PassbandFrequency = 3500;
[num,den] = tf(LPF);

% fvtool(num,den);


Signal_time = real(ifft(Signal));
Signal_time = conv(Signal_time, num);
Signal_time = Signal_time / (max(Signal_time));
Signal_time = circshift(Signal_time, round(length(Signal_time) * 0.3));
Signal_freq = fft(Signal_time);


figure(1);
subplot(1, 2, 2);
semilogx(f, 10 * log10(abs(Signal_freq(1 : round(length(Signal_freq) / 2)))));
title("sweep band signal frequency domain");
xlabel('[Hz]');

subplot(1, 2, 1);
plot(Signal_time);
ylim([-1.5 1.5]);
title("sweep band signal time domain");
xlabel('step');
save('Gaussian_band.mat', 'Signal_time');


