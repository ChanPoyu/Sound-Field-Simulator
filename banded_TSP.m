clear all
close all

%% set all parameters
% l = 13;
fs = 48000;
N = 10000;
m = N / 4;

%% swept sine signal in frequency domain

S = zeros(1, N);

for n = 0 : 1 : N/2
    
    S(n + 1) = exp((-1) * 1i * 4 * m * pi * n^2 / N^2);
    
    
%         omega = n / (N / 2) * (fs / 2) * 2 * pi ;
%         weight = 1 / (1i * omega ^ 1.5);
%         S(n+1)  = S(n+1) * weight;
   
    
end
S(1) = 0;
for n = N/2 + 1 : 1 : N-1
    S(n + 1) = conj(S(N - n + 1));
end



%% inverse fft
fupper = 4000;
flower = 100;
fs = 48000;

% band_bottom = floor(flower / (fs/2) * 2500);
% band_top = ceil(fupper / (fs/2) * 2500);
% S(1:band_bottom) = 0;
% S(N-band_bottom+1:end) = 0;
% S(band_top:N-band_top+1) = 0;

% for n = 11 : 1 : 838
%     S(n) = (n - 10) / 828;
% end
% for n = 839 : 1 : 1665
%     S(n) = (1666 - n) / 828;
% end
% for n = 3336 : 1 : 4163
%     S(n) = (n - 3335) / 828;
% end
% for n = 4164 : 1 : 4990
%     S(n) = (4991 - n) / 828;
% end

TSP = ifft(S);
TSP = circshift(TSP, N/4);
TSP = TSP / max(TSP);

save('TSP.mat','TSP');

%% bandpass filter
% f = [0 100/48000 2000/48000 1];    % point of frequency window
% m = [0 1 1 0];                     % amplitude
% band_pass_filter = fir2(30, f, m); % generate fir filter
% [h1, w] = freqz(band_pass_filter, 1);
% 
% figure(10);
% plot(f, m, w/pi, abs(h1));

%% lowpassfilter
LPF = dsp.LowpassFilter;
LPF.SampleRate = 48000;
LPF.StopbandFrequency = fupper * 1.25;
LPF.PassbandFrequency = fupper;
[num,den] = tf(LPF);
fvtool(num,den);

%% convolution
TSP_filtered = conv(TSP, num);
TSP_filtered = TSP_filtered - TSP_filtered(1);
TSP_filtered = TSP_filtered / max(TSP_filtered);
Y_filtered = fft2(TSP);

save('TSP_filtered.mat', 'TSP_filtered');

%% visualize results
figure(1);
% subplot(2,2,1);
% plot(1:N, real(TSP), 1:N, imag(TSP));
% title('TSP [time domain]');
% 
% 
% subplot(2,2,2);
% semilogx(1 : fs/N :fs/2 , 20 * log10(abs(S(1 : length(S) / 2))));
% title('TSP [frequency domain]');


subplot(1,2,1);
plot(real(TSP_filtered));
title('TSP Filtered [time domain]');
axis tight;


subplot(1,2,2);
semilogx(1 : fs/N :fs/2, 20 * log10(abs(Y_filtered(1 : length(Y_filtered) / 2))));
ylim([0 50]);
title('TSP Filtered [frrequency domain]');
axis tight;

audiowrite('TSP.wav', TSP, 48000);
audiowrite('TSP_filtered.wav', TSP_filtered, 48000);



%% swept sine invers

S_inverse = zeros(1, N);

for n = 0 : 1 : N/2
  S_inverse(n + 1) = exp(1i * 4 * m * pi * n^2 / N^2);
end

for n = N/2 + 1 : 1 : N-1
    S_inverse(n + 1) = conj(S_inverse(N - n + 1));
end

% cut off without flower ~ fupper

TSP_inverse = ifft(S_inverse);
TSP_inverse = circshift(TSP_inverse, N * 3 / 4);
TSP_inverse = TSP_inverse / max(TSP_inverse);
TSP_inverse_filtered = conv(TSP_inverse, num);

save('TSP_inverse.mat','TSP_inverse');
save('TSP_inverse_filtered.mat','TSP_inverse_filtered');

figure(50);
plot(1:N, real(TSP_inverse), 1:N, imag(TSP_inverse));



inpulse_lowpass = conv(TSP, TSP_inverse);
figure(52);
plot(real(inpulse_lowpass));
title('inpulse lowpassed');

