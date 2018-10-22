close all
clear all

%% parameter to design a TSP
fs = 48000;
f_upper = 4000;
f_lower = 100;

N = 5000;
J = N * 1 / 2;

f = linspace (0, fs / 2, floor(N /2) + 1);

cut_freq = 100;  % -3db from cut frequency

%% f characteristic of MAGNITUDE in [db]

magnitude = ones(1 , length(f));

% for n = 1 : 1 : length(magnitude)
%     magnitude(n) = 10 * log10(1 / (1 + (f(n) / cut_freq) ^ 2)) - 3;
% end


envelope = ones(1, length(f));

% for n = 1 : 1 : length(envelope)
%    envelope(n) = 10 * log10(1 / (1 + (f(n) ^ 1.5 / cut_freq) ^ 2)) - 3;
% end
;

E = zeros(1, length(magnitude));

for n = 1 : 1 : length(E)
    E(n) = (10 ^ (magnitude(n) / 10)) / (10 ^ (envelope(n) / 10));
end

% group delay
gd = cumsum(E);
gd = gd * (J / gd(end));

% phase
ph = cumsum(gd) * 2 * pi / N;
ph = ph * round(ph(round(N / 2)) / pi) * pi / ph(round(N / 2));

% amplitude
A = zeros(1, length(E));
for n = 1 : 1 : length(A)
    A(n) = sqrt(E(n) * (J * N) / 2 / (2 * sum(E) - E(1) - E(round(N / 2))));
    A(n) = A(n) * 10 ^(envelope(n) / 20);
end

% IFFT
SS = zeros(1, N);

for n = 1 : 1 : N / 2 + 1
   SS(n) = A(n) * exp( -1i * ph(n)); 
end

for n = N / 2 + 2 : 1 : length(SS)
    SS(n) = conj(SS(N - n + 2));
end


figure(100);
plot(real(SS));
hold on;
plot(imag(SS));
hold off;



ss = real(ifft(SS));
ss = circshift(ss, round((N-J) / 2));
ss = ss / max(ss);

%% arbitrary inverse

ss_inverse = real(ifft(1 ./ fft(ss)));
% ss_inverse = circshift(ss_inverse, round((N-J) / 2));
ss_inverse = ss_inverse / max(ss_inverse);

SS_inverse = fft(ss_inverse);


figure(2);
subplot(2,1,1);
semilogx(f, abs(SS_inverse(1: length(f))));

subplot(2,1,2)
plot(real(ss_inverse));


audiowrite('arbitrary_TSP_inverse.wav', ss_inverse, 48000);
save('arbitrary_TSP_inverse.mat', 'ss_inverse');






%% LPF

LPF = dsp.LowpassFilter;
LPF.SampleRate = 48000;
LPF.StopbandFrequency = f_upper * 1.25;
LPF.PassbandFrequency = f_upper;
[num,den] = tf(LPF);
% fvtool(num,den);

% ss_db = 10 * log10(real(ss) .^ 2);
% figure(5);
% plot(ss_db);
% hold on;


% 
% ss_db = 10 * log10(real(ss) .^ 2);
% plot(ss_db);


%% 1 / omega filter

FIR = sqrt(1 ./ (1 + (f / cut_freq) .^ (2 * 1.5)));

wind = blackmanharris(1027)';
fir = fir2(1025, f / (fs / 2), FIR);

figure(11);
subplot(2, 1, 1)
plot(fir);

subplot(2, 1, 2);
semilogx(10 * log10(abs(fft(fir, 48000))));


ss = conv(num, ss);

pulse = conv(ss, ss_inverse);
pulse_freq = fft(pulse);

ss = conv(fir, ss);




figure(1);

subplot(2, 2, 1);
semilogx(f, magnitude);
title("mag. sweep");
xlabel("Hz");
ylabel("db");
axis tight;


subplot(2, 2, 2);
semilogx(f, envelope);
title("desire envelope");
xlabel("Hz");
ylabel("db");
axis tight;

subplot(2, 2, 3);
semilogx(f, gd / fs);
title("delay sweep");
xlabel("Hz");
ylabel("sec");
axis tight;

subplot(2, 2, 4);
% plot(0 : 1 / fs :(length(ss) - 1) / fs, real(ss));
plot(real(ss));
title("sweep");
xlabel("sec");
ylabel("V");
axis tight;

figure(3);
% subplot(2, 1, 1);
plot(pulse);
axis tight;
% 
% subplot(2, 1 ,2);
% plot(abs(pulse_freq(1 : round(end / 2))));


audiowrite('arbitrary_TSP.wav', ss, 48000);
save('arbitrary_TSP.mat', 'ss');




