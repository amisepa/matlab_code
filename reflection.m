%% Implementation of the reflection method by Mike X Cohen to avoid undesired 
% edge effects when performing time-frequency analysis on EEG signals. 
% 
% Cedric Cannard, 2023

%% create a signal

N = 500;
hz = linspace(0,1,N);

% create a gaussian (see time domain denoising video)
gx = exp( -(4*log(2)*(hz-.1)/.1).^2 )*N/2;  

% generate random phase values between 0 and 2*pi are generated (rand part), 
% insert them into Euler's formula to give complex unit-length random factors, 
% point-wise multiply them by the gaussian (gx)
% and take the inverse fourier transform
% this creates a signal that has a lot of random values within some
% frequency range
signal = real(ifft( gx.*exp(1i*rand(1,N)*2*pi) )) + randn(1,N); 

% Plot it and its power
figure(1), clf
subplot(311)
plot(1:N,signal,'k')
set(gca,'xlim',[0 N+1])
title('Raw simulated signal')
xlabel('Time points (a.u.)')

subplot(324)
plot(hz,abs(fft(signal)).^2,'k','markerfacecolor','w')
set(gca,'xlim',[0 .5])
xlabel('Frequency (norm.)'), ylabel('Energy')
title('Frequency-domain signal representation')

%% Apply a low-pass zero-phase filter

% generate filter kernel
order = 150;
fkern = fir1(order,.6,'low'); %.6 is in term sof the Nyquist frequency C9orrespond to .3 on the frequency axis of the plot)
% figure; plot(fkern)

% zero-phase-shift filter
fsignal = filter(fkern,1,signal);               % forward
fsignal = filter(fkern,1,fsignal(end:-1:1));    % reverse
fsignal = fsignal(end:-1:1);                    % flip forward

% plot the original signal and filtered version
subplot(323), hold on
plot(1:N,signal,'k')
plot(1:N,fsignal,'m')
set(gca,'xlim',[0 N+1])
xlabel('Time (a.u.)')
title('Time domain')
legend('original','filtered (no reflection)')
xlabel('Frequency (norm.)'), ylabel('Energy')

% power spectra
subplot(324), hold on
plot(hz,abs(fft(fsignal)).^2,'m')
title('Frequency domain')
legend('Original','filtered (no reflection)')
set(gca,'xlim',[0 .5])

% Here, we can see edge effect in the time domain (flat line) and the power 
% spectrum shoulf be the same curve as before but without power above the cutoff
% at .3 on x-axis, but it is deformed

%% apply the reflection method by filter order

% reflect the signal by adding a backward-version of the signal before and
% after. 
reflectsig = [ signal(order:-1:1) signal signal(end:-1:end-order+1) ]; % here we use the length of the filter kernel as the length of the reflection
% refectsig = [ signal(end:-1:1) signal signal(end:-1:1) ]; % here we use the length of the whole signal 

% zero-phase-shift filter the reflected signal
reflectsig = filter(fkern,1,reflectsig);              % forward
reflectsig = filter(fkern,1,reflectsig(end:-1:1));    % reverse
reflectsig = reflectsig(end:-1:1);                    % flip forward

% chop off the reflected parts (now it has a length of 500 like the
% original signal)
fsignal = reflectsig(order+1:end-order);

% plot
subplot(325), hold on
plot(1:N,signal,'k')
plot(1:N,fsignal,'m')
set(gca,'xlim',[0 N+1])
xlabel('Time (a.u.)');
title('Time domain')
subplot(326), hold on
plot(hz,abs(fft(signal)).^2,'k')
plot(hz,abs(fft(fsignal)).^2,'m')
title('Frequency domain')
legend('Original','filtered (with reflection)')
set(gca,'xlim',[0 .5])
xlabel('Frequency (norm.)'), ylabel('Energy')


%% try again with filtfilt (requires signal processing toolbox)

fsignal1 = filtfilt(fkern,1,signal);

figure
plot(1:N,fsignal,1:N,fsignal1)
%  --> see only one line because filtfilt applies the reflection, but minor 
% differences if we zoom in, because of some additional padding

corr(fsignal(:),fsignal1(:))
