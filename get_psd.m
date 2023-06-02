%% Compute power spectral density (PSD) or power for each EEG channel using
% MATLAB's pwelch method.  Defaults are Hamming taper on a 2-s window with
% 50% overlap, outputting the power spectral density (PSD).
% 
% Usage:
% [psd, freqs] = get_psd(eeg_data,winSize,taperM,overlap,nfft,Fs,freqRange,type,vis);
% [psd, freqs] = get_psd(EEG.data,EEG.srate*2,'hamming',50,[],EEG.srate,[1 100],'psd',1);
% 
% - eeg_data with channels in 1st dimension and data in 2nd dimension (default = EEG.data)
% - window size in frames (default = 2 s window).
% - taper method: hamming (default), hann, blackman, rectwin.
% - overlap in percent (default = 50)
% - Fs: sample rate in Hz (default = EEG.srate)
% - freqRange is frequecnies of interest to compute (default = 1:100)
% - type: returns power spectral density ('psd'; default) or returns
%           'power' (scales each estimate of the PSD by the equivalent noise 
%           bandwidth of the window (in hertz): i.e. power estimate at each frequency).
% 
% Cedric Cannard, 2021

function [pwr, pwr_norm, f] = get_psd(eegData,winSize,taperM,overlap,nfft,Fs,fRange,type,vis)

% Error if no sampling rate provided
if ~exist('Fs', 'var') || isempty(Fs)
    errordlg('You need to provide the sampling rate Fs to use this function.'); return;
end

% Window size
if ~exist('winSize', 'var') || isempty(winSize) 
    disp('Window size not provided: 2-s windows');
    winSize = Fs*3;
end

% Taper
if ~exist('taperM', 'var') || isempty(taperM) 
    taperM = 'hamming'; %hamming (default); hann; blackman; rectwin
end
fh = str2func(taperM);

% Overlap
if ~exist('overlap', 'var') || isempty(overlap)
    overlap = 50;
end
overlap = winSize/(100/overlap); % convert overlap to samples

% Frequency range default
if ~exist('fRange', 'var') || isempty(fRange)
    nyquist = Fs/2;
    fRange = [1/nyquist nyquist]; 
end

% Power type default
if ~exist('type', 'var')
    type = 'psd';    
end

% nfft
if ~exist('nfft', 'var') 
    nfft = Fs*2;
end

if vis
    figure('color','w'); hold on
end

% Power spectral density (PSD)
for iChan = 1:size(eegData,1)
    [pwr(iChan,:), f] = pwelch(eegData(iChan,:),fh(winSize),overlap,nfft,Fs,type);
    if vis
        plot(f,log10(pwr(iChan,:))); 
    end
end

% Calculate frequency resolution
% exp_tlen = nextpow2(tlen);
% fres = Fs/2.^exp_tlen;

% Truncate PSD to frequency range of interest (ignore freq 0)
% freq = dsearchn(f,fRange(1)):dsearchn(f, fRange(2));
% f = f(freq(2:end))';
freq = f >= fRange(1) & f <= fRange(2);
f = f(freq);
pwr = pwr(:,freq(2:end));     

% Normalize to deciBels (dB)
pwr_norm = 10*log10(pwr);

