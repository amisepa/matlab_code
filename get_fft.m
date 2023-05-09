%% get_coherence
% Compute power spectrum on each sliding window using fft.
%  
% Inputs:
%   - EEG data (continuous)
%   - sRate: data sampling rate in Hz
%   - wSize: window size in s (default = 2)
%   - freqRange: frequency range to output (e.g [1 50]). 
%               default = 0 to Nyquist limit.
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

function [p_mu, p, f] = get_fft(data, fs, winSize, freqRange, norm)

dt = 1/fs;   %discrete sampling in s (0.002 = 2 ms)
df = 1/winSize;   %frequency resolution
fNQ = 1/dt/2;   %Nyquist frequency
f = df:df:fNQ;  %whole frequency vector availabe in data

% error if no sampling rate provided
if ~exist('fs', 'var') || isempty(fs)
    errordlg('You need to provide the sampling rate Fs to use this function.'); return;
end

% Window size
if ~exist('winSize', 'var') || isempty(winSize) 
    disp('Window size not provided: 2-s windows');
    winSize = fs*2;
end

% Frequency range    
if ~exist('freqRange', 'var') || isempty(freqRange)
    disp('Computing for all frequencies up to Nyquist limit');
    freqRange = f;
end

% Normalize to dB 
if ~exist('norm', 'var') || isempty(norm)
    norm = false;
end

% Compute power spectra for each sliding window
nWind = floor((size(data,2)/fs)/winSize);
epoch = 1;
p = zeros(size(data,1),1,fs*winSize); %preallocate memory
for iChan = 1:size(data,1)
    for iWind = 1:winSize*fs:nWind*fs-1%winSize*fs
        p(iChan,epoch,:) = 2*dt^2/winSize * fft(data(iChan,iWind:iWind+winSize*fs-1)).*conj(fft(data(iChan,iWind:iWind+winSize*fs-1)));
        epoch = epoch+1;
    end
end

% Keep positive frequencies
p = p(:,:,1:size(p,3)/2);

% Keep only frequencies of interest
if size(freqRange,2) == 2
    idx = and( f>=freqRange(1), f<=freqRange(2) );
    p = p(:,:,idx);
    % f = freqRange(1):1/winSize:freqRange(end);
    f = f(:,idx);
elseif size(freqRange,2) == 1
    idx = f==freqRange; %if only one frequency was selected
    p = p(:,:,idx);
    f = f(:,idx);
end

% Normalize to dB
if norm
    p = 10*log10(p);
end

% Average across epochs
p_mu = squeeze(mean(p,2))';

