%% Compute power spectral density (PSD) or power for each EEG channel using
% MATLAB's pwelch method.  Defaults are Hamming taper with 50% overlap.
% 
% Usage:
% [psd, freqs] = get_psd(EEG.data, frames, taperM, overlap, nfft, Fs, freqRange);
% [psd, freqs] = get_psd(EEG.data, EEG.srate, 'hamming', 50, [], EEG.srate, 1:50);
% 
% - eegData with channels as 1st dimension and data as 2nd dimension (default = EEG.data)
% - number of frames per window (default = sample rate for 1 s window)
% - taper method: hamming, hann, blackman, rectwin (default is 'hamming').
% - overlap in percent (default = 50)
% - Fs is sample rate in Hz (default = EEG.srate)
% - freqRange is frequecnies of interest to compute (default = 1:100)
% - type: returns power spectral density ('psd'; default) or returns
%           'power' (scales each estimate of the PSD by the equivalent noise 
%           bandwidth of the window (in hertz): i.e. power estimate at each frequency).
% 
% Cedric Cannard, 2021

function [pxx, f] = get_psd(eegData, frames, taperM, overlap, nfft, Fs, fRange)

%error if no sampling rate provided
if ~exist('Fs', 'var') || isempty(Fs)
    errordlg('You need to provide the sampling rate Fs to use this function.'); return;
end

%Frames
if ~exist('frames', 'var') || isempty(frames) 
    if ~exist('Fs', 'var') || ~isempty(frames) 
        disp('Frames not provided: taking sample rate for 1 s windows');
        frames = Fs;
    end
end

%Tapering
if ~exist('taperM', 'var') || isempty(frames) 
    taperM = 'hamming'; %hann; blackman; 
end
% if strcmp(taperM, 'hamming')
%     window = hamming(frames);
% elseif strcmp(taperM, 'hann')
%     window = hann(frames);
% elseif strcmp(taperM, 'blackman')
%     window = blackman(frames);
% elseif strcmp(taperM, 'rectwin')
%     window = rectwin(frames);
% end
fh = str2func(taperM);

% Overlap
if ~exist('overlap', 'var') || isempty(overlap)
    overlap = 50;
end
overlap = frames/(100/overlap); %convert overlap in samples for pwelch

% Frequency range default
% if ~exist('fRange', 'var') || isempty(fRange)
%     fRange = [1 Fs/2];     %1 Hz for the min and Nyquist (srate/2) for the max. 
% end

% power type default
% if ~exist('type', 'var') || isempty(type)
%     type = 'psd';     %1 Hz for the min and Nyquist (srate/2) for the max. 
% end

%nfft
if ~exist('nfft', 'var') 
    nfft = [];
end

% %Detrend data
% if size(eegData,3) > 1
%     for iTrial = 1:size(eegData,3)
%         eegData(:,:,iTrial) = detrend(eegData(:,:,iTrial)')';
%     end
% else
%     eegData = detrend(eegData')';
% end
    
%Power spectral density (PSD)
% disp('computing power spectral density (PSD) for each channel...');
for iChan = 1:size(eegData,1)
%     [pxx(:,iChan), f] = pwelch(eegData(iChan,:), window, overlap, 1:fRange(end), Fs, type);
    [pxx(:,iChan), f] = pwelch(eegData(iChan,:), fh(frames), overlap, nfft, Fs);
%     [pxx(:,iChan), f, c] = pwelch(eegData(iChan,:), window, overlap, [], Fs, type);
%     [pxx(:,iChan), f] = pwelch(data(iChan,:), window, overlap, nfft, Fs);
end

%Calculate frequency resolution
% exp_tlen = nextpow2(tlen);
% fres = Fs/2.^exp_tlen;

frex = dsearchn(f, fRange(1)):dsearchn(f, fRange(2));
f = f(frex);
pxx = pxx(frex,:)';     % truncate PSD to frex range
% pxx = 10*log10(pxx);    %normalize to dB

% freqs = round(f)>=fRange(1) & round(f)<=fRange(end);  %get correct frequencies
% power = 10*log10(pxx(f,:))';      %Convert to dB and keep only freqs of interest
% psd = log(pxx(freqs,:))';         %Convert to natural log (better for differenciation)
% c = 10*log10(c(freqs,:))';        %confidence bounds
% freqs(freqs == 0) = [];
% freqs = find(freqs);

% psd = pxx(freqs,:);
% figure; plot(freqs, psd);hold on; plot(freqs,log(psd));hold on; plot(freqs, 10*log10(psd));
% legend('psd', 'log', 'dB', 'power');

% Frequencies in Welch's estimate where the lower confidence bound exceeds 
% the upper confidence bound for surrounding frequencies clearly indicate 
% significant oscillations in the time series.Default is 95% confidence bounds
% color1 = [0, 0.4470, 0.7410];          
% figure; set(gcf,'Color','w');
% plot(freqs', psd,'LineWidth',2,'Color',color1); hold on;
% fillhandle = fill([freqs' fliplr(freqs')],[c(1,:) fliplr(c(2,:))], color1);
% set(fillhandle,'EdgeColor', color1,'FaceAlpha',0.3,'EdgeAlpha',0.8);

end
