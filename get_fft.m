%% get_coherence
% Compute power spectrum on each sliding window using fft.
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
%   [p_mu, p, f] = get_fft(data,fs, 2, 8:13)
%   figure; plot(f,p_mu); xlabel('Frequencies (Hz)'); ylabel('Power spectrum');
%   or 
%   p_mu = get_fft(data, EEG.srate, 2, 10)
% 
% Cedric Cannard, January 2022

function [p_mu, p, f] = get_fft(data, fs, winSize)

dt = 1/fs;   %discrete sampling in s (0.002 = 2 ms)
df = 1/winSize;   %frequency resolution
fNQ = 1/dt/2;   %Nyquist frequency
f = df:df:fNQ;  %whole frequency vector availabe in data

%error if no sampling rate provided
if ~exist('fs', 'var') || isempty(fs)
    errordlg('You need to provide the sampling rate Fs to use this function.'); return;
end

%Window size
if ~exist('winSize', 'var') || isempty(winSize) 
    disp('Window size not provided: 2-s windows');
    winSize = fs*2;
end

% %check frequency vector    
% if ~exist('freqs', 'var') || isempty(winSize) 
%     disp('Computing for all frequencies up to Nyquist limit');
%     freqs = f;
% end

%Initiate variables
nWind = floor((size(data,2)/fs)/winSize);
p = zeros(size(data,1),1,fs*winSize); %preallocate memory
epoch = 1;

%Compute power spectra for each sliding window
for iChan = 1:size(data,1)
    for iWind = 1:winSize*fs:nWind*fs-1%winSize*fs
        p(iChan,epoch,:) = 2*dt^2/winSize * fft(data(iChan,iWind:iWind+winSize*fs-1)).*conj(fft(data(iChan,iWind:iWind+winSize*fs-1)));
        epoch = epoch+1;
    end
end

%Keep positive frequencies
p = p(:,:,1:size(p,3)/2);

% % Keep only frequencies of interest
% if size(freqs,2) > 1
%     f = find(f==freqs(1)):find(f==freqs(end));
%     p = p(:,f);
%     f = freqs(1):1/winSize:freqs(end);
% else
%     f = find(f==freqs); %if only one frequency was selected
%     p = p(:,f);
% end

%Average across epochs
p_mu = mean(p,2);
if size(data,1)>1
    p_mu = squeeze(p_mu);
end

