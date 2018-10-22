close all;
clear all;

dt = 0 : 1 / 48000 : 10000 / 48000 ;

P = sin(2 * pi * 500 * dt) .^ 2;

f = linspace(0 ,  24000, 5000);

Y = fft(P);
figure(1);

subplot(1,2,1);
plot(P(1 : 1000));
title('sine ^2 in time domain');
xlabel('step');
axis tight;

subplot(1,2,2);
semilogx(f, 10 * log10(abs(Y(1 : (end-1)/2))));
title('Gaussian Pulse in frequency domain');
xlabel('[Hz]');

axis tight;