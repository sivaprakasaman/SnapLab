clear all;
close all;
%% Parameters:
Fs = 48828.125;
sec = 1.3; %stim length in seconds
f_carrier = 4000; %carrier frequency
f_mod = 100; %modulation freq
window =[0.45,0.55]; %window for plotting
duty1 = 50; %duty cycle (max 1 = 100%)
duty2 = 25; %duty cycle 2

%% SAM Tone:
t = 0:1/Fs:sec; %samples
t = t(:);
%4kHz Carrier Frequency:
% = sin(f_carrier*2*pi*t);
s_carrier = sin(f_carrier*2*pi*t);

%120Hz Modulation Frequency:
s_mod = 0.5 + 0.5.*sin(f_mod*2*pi*t);

%SAM tone Generation:
SAM = s_mod.*s_carrier;
% plot((window(1)*Fs:window(2)*Fs)/Fs,SAM(window(1)*Fs:window(2)*Fs))%first half-second of SAM
%% Square Tone
sq = square(2*pi*f_mod*t,duty1);
sq = sq.*(sq>=0);
sq_50 = sq.*s_carrier;
sq_1 = sq_50;

sq = square(2*pi*f_mod*t,duty2);
sq = sq.*(sq>=0);
sq_25 = sq.*s_carrier;
sq_2 = sq_25;

% figure;
% plot((window(1)*Fs:window(2)*Fs)/Fs,sq_120(window(1)*Fs:window(2)*Fs))%first half-second of SAM
%% 10 H Complex

H(1,:) = sin(f_mod*2*pi*t).*s_carrier;
H(2,:) = sin(2*f_mod*2*pi*t).*s_carrier;
H(3,:) = sin(3*f_mod*2*pi*t).*s_carrier;
H(4,:) = sin(4*f_mod*2*pi*t).*s_carrier;
H(5,:) = sin(5*f_mod*2*pi*t).*s_carrier;
H(6,:) = sin(6*f_mod*2*pi*t).*s_carrier;
H(7,:) = sin(7*f_mod*2*pi*t).*s_carrier;
H(8,:) = sin(8*f_mod*2*pi*t).*s_carrier;
H(9,:) = sin(9*f_mod*2*pi*t).*s_carrier;
H(10,:) = sin(10*f_mod*2*pi*t).*s_carrier;


H10 = sum(H)/max(sum(H));

%plot(sum(H)./max(sum(H)));

%% Plotting

subplot(4,1,1)
plot((window(1)*Fs:window(2)*Fs)/Fs,SAM(window(1)*Fs:window(2)*Fs))%first half-second of SAM
title('SAM');

subplot(4,1,2)
plot((window(1)*Fs:window(2)*Fs)/Fs,sq_1(window(1)*Fs:window(2)*Fs))%first half-second of SAM
title('50% Duty Cycle')

subplot(4,1,3)
plot((window(1)*Fs:window(2)*Fs)/Fs,sq_2(window(1)*Fs:window(2)*Fs))%first half-second of SAM
title('25% Duty Cycle');
xlabel('Time (s)')

subplot(4,1,4)
plot((window(1)*Fs:window(2)*Fs)/Fs,H10(window(1)*Fs:window(2)*Fs))%first half-second of SAM
title('10 H complex');
xlabel('Time (s)')

%% Play stimulus
player = audioplayer([sq_25,sq_25],Fs);
player.play();

%% Save Files
save('sq_50','sq_50');
save('sq_25','sq_25');
save('SAM','SAM');

%% FFT?
sig2 = sq_25;

fft_sig2 = fft(sig2);

T = 1/Fs;             % Sampling period       
L = length(sig2);     % Length of signal
t = (0:L-1)*T;        % Time vector

P2 = abs(fft_sig2/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure;
subplot(311)
plot(f,P1) 
ylabel("Amplitude")
title("Spectrum of SQ_25 signal")


sig2 = sq_50;

fft_sig2 = fft(sig2);

T = 1/Fs;             % Sampling period       
L = length(sig2);     % Length of signal
t = (0:L-1)*T;        % Time vector

P2 = abs(fft_sig2/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
subplot(312)
plot(f,P1) 
xlabel("Frequency")
ylabel("Amplitude")
title("Spectrum of SQ50 signal")

sig2 = SAM;

fft_sig2 = fft(sig2);

T = 1/Fs;             % Sampling period       
L = length(sig2);     % Length of signal
t = (0:L-1)*T;        % Time vector

P2 = abs(fft_sig2/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
subplot(313)
plot(f,P1) 
xlabel("Frequency")
ylabel("Amplitude")
title("Spectrum of SAM signal")

