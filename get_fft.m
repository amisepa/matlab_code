%% get_coherence
% Compute power-spectrum on each sliding window using fft.
%  
% Inputs:
%   - EEG data (continuous)
%   - sRate: data sampling rate in Hz
%   - wSize: window size in s (default = 2)
%   - fInt: vector of frequencies of interest for coherence output; default = up to Nyquist frequency)
% 
% Outputs:
%   - p_mu: power spectrum averaged across windows
%   - p: power spectrum for each window
%   - f: frequency vector
% 
% Example: 
%   [p_mu, p, f] = get_fft(data, EEG.srate, 2, 8:13)
%   figure; plot(f,p_mu); xlabel('Frequencies (Hz)'); ylabel('Power spectrum');
%   or 
%   p_mu = get_fft(data, EEG.srate, 2, 10)
% 
% Cedric Cannard, January 2022

function [p_mu, p, f] = get_fft(data, sRate, wSize, fInt)

dt = 1/sRate;   %discrete sampling in s (0.002 = 2 ms)
df = 1/wSize;   %frequency resolution
fNQ = 1/dt/2;   %Nyquist frequency
f = df:df:fNQ;  %whole frequency vector availabe in data

%check frequency vector    
if fInt(1)<1
    error('Minimum frequency should be 1 Hz.');
end
if fInt(end)> f(end)
    warning('Upper frequency bound selected is above the Nyquist limit. Changing to Nyquist frequency.');
    fInt = fInt(1):f(end);
end

%Initiate variables
nWind = floor((size(data,2)/sRate)/wSize);
p = zeros(1,sRate*wSize);
epoch = 1;

%Compute power spectra and cross-spectrum for each sliding window
for iWind = 1:wSize*sRate:nWind*sRate-wSize*sRate
    p(epoch,:) = 2*dt^2/wSize * fft(data(iWind:iWind+wSize*sRate-1)).*conj(fft(data(iWind:iWind+wSize*sRate-1)));
    p_sd(epoch,:) = std(p(epoch,:));
    epoch = epoch+1;
end

%Keep positive frequencies
p = p(:,1:size(p,2)/2);

% Keep only frequencies of interest
if size(fInt,2) > 1
    f = find(f==fInt(1)):find(f==fInt(end));
    p = p(:,f);
    f = fInt(1):1/wSize:fInt(end);
else
    f = find(f==fInt); %if only one frequency was selected
    p = p(:,f);
end

%Average across epochs
p_mu = mean(p,1);
