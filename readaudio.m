clear all
close all

%% read audio

f_upper = 4000;

[clap, fs] = audioread('pulse.wav');



LPF = dsp.LowpassFilter;
LPF.SampleRate = fs;
LPF.StopbandFrequency = f_upper * 1.25;
LPF.PassbandFrequency = f_upper;
[num,den] = tf(LPF);
% fvtool(num,den);

lowpass_clap = conv(num, clap);
lowpass_clap = lowpass_clap(63000: 75000);
CLAP = fft(lowpass_clap);
save('clap_lowpass.mat', 'lowpass_clap');


figure(1);
subplot(2,1,1);
plot(clap);

subplot(2,1,2);
plot(lowpass_clap);

figure(2);
semilogx(10 * log10(abs(CLAP(1 : round(length(CLAP) /2)))));