%% Compute power spectrum for each channel available in data.
% 
% EXAMPLE: 
%   [pwr, pwr_norm, psd, psd_norm, f, pwr_seg, psd_seg] = compute_pwr(data, fs, overlap, fRange, winLength, coherence,vis)
%   [pwr, pwr_norm, psd, psd_norm, f, pwr_seg, psd_seg] = compute_pwr(EEG.data, EEG.srate, .5, [1 100], 2, true, true)
% 
% Copyright (C) - Cedric Cannard, 2024

function [pwr, pwr_norm, psd, psd_norm, f, psd_seg, t_seg, coh_seg, pcoh_seg, gpdc_seg, coh_f] = compute_pwr(data, fs, overlap, fRange, winLength, coherence, vis)

coh_seg = [];
pcoh_seg = [];
gpdc_seg = [];
coh_f = [];

% check data dimensions, ensure channels x samples
if size(data,1) > size(data,2)
    data = data';
end

% Sampling rate
if ~exist('fs', 'var') || isempty(fs)
    errordlg('You need to provide the sampling rate Fs to use this function.'); return;
end

% Epoched vs continuous
if length(size(data)) == 2
    epoched = false;
    disp('Continuous data detected.')
else
    epoched = true;
    disp('Epoched data detected. Converting them to continuous.')
    data = data(:,:);
end

nChan = size(data,1);
nSamples = size(data,2);

% Window size
if ~exist('winLength','var')
    winLength = 2;   % 2-s window by default
end
winSize = fs * winLength;  % in samples

% Taper method: hamming (default)
win = hamming(winSize);  % Create the Hamming window

% Overlap
if ~exist('overlap', 'var') || isempty(overlap)
    overlap = .5;
end
noverlap = floor(overlap * winSize);  % Convert overlap ratio to samples

% Frequency range default
if ~exist('fRange', 'var') || isempty(fRange)
    nyquist = fs / 2;
    fRange = [1 nyquist]; 
end

% nfft
if ~exist('nfft', 'var') 
    nfft = winSize*2;  % Number of FFT points
end

% Su
if coherence
    Su = eye(nChan, nChan);
end

% Segment the data with 50% overlap
step = winSize - noverlap;  % Step size
segIdx = 1:step:(nSamples - winSize + 1);   % Starting indices for segments
nSegments = length(segIdx);                 % Number of segments
seg_end_idx = segIdx + winSize - 1;         % End index of each window
seg_end_time = seg_end_idx / fs;            % Convert to seconds

% Perform FFT and calculate PSD for each windowed segment across all channels
psd_seg = nan(nChan, floor(nfft / 2) + 1, nSegments);
pwr_seg = nan(nChan, floor(nfft / 2) + 1, nSegments);
coh_seg = nan(nChan, nChan, nfft, nSegments);
pcoh_seg = nan(nChan, nChan, nfft, nSegments);
gpdc_seg = nan(nChan, nChan, nfft, nSegments);
if ~coherence
    progressbar('Computing spectral power over sliding windows')
else
    progressbar('Computing spectral power  + coherence over sliding windows')
end
for iSeg = 1:nSegments
    idx_start = segIdx(iSeg);
    idx_end = idx_start + winSize - 1;
    t_seg(iSeg) = seg_end_time(iSeg); % end time in seconds

    segment = data(:, idx_start:idx_end);  % [nChan x winSize]
    segment = segment .* win';             % Apply window tapering

    Y = fft(segment, nfft, 2);             % FFT across time dimension

    psd_seg(:, :, iSeg) = (1 / (fs * winSize)) * abs(Y(:, 1:floor(nfft/2) + 1)).^2; % power spectral density
    pwr_seg(:, :, iSeg) = abs(Y(:, 1:floor(nfft/2) + 1)).^2;

    if coherence
        [dc, dtf, pdc, gpdc, ~, coh, pcoh, ~, ~, ~, ~, coh_f] = fdMVAR_5order(segment, Su, nfft, fs);
        coh_seg(:,:,:,iSeg) = abs(coh);     % Coherence between two signals at each frequency
        pcoh_seg(:,:,:,iSeg) = abs(pcoh);   % Partial Coherence (after removing the influence of all other signals, isolating unique shared variance)
        gpdc_seg(:,:,:,iSeg) = abs(gpdc);   % Direct influence from one signal to another, normalized by total output from the source, accounting for noise levels in each signal.
    end

    progressbar(iSeg / nSegments)
end

% Average PSD and power across segments
psd = trimmean(psd_seg, 20, 3);  % average across segments
pwr = trimmean(pwr_seg, 20, 3);  % average across segments

% Frequency vector
f = fs * (0:(floor(nfft / 2))) / nfft;

% Keep only frequencies of interest (ignore freq 0)
freq_idx = f >= fRange(1) & f <= fRange(2);
f = f(freq_idx);
psd = psd(:, freq_idx);     
pwr = pwr(:, freq_idx);
psd_seg = psd_seg(:,freq_idx,:);
% pwr_seg = pwr_seg(:, mask);

% Same for coherence data
if coherence
    freq_idx = coh_f >= fRange(1) & coh_f <= fRange(2);
    if coh_f(freq_idx(1)) == 0, freq_idx(1) = []; end
    coh_f = coh_f(freq_idx);
    coh_seg = coh_seg(:,:,freq_idx,:);
    pcoh_seg = pcoh_seg(:,:,freq_idx,:);
    gpdc_seg = gpdc_seg(:,:,freq_idx,:);
end

% Normalize power to decibels (dB)
psd_norm = 10 * log10(psd);
pwr_norm = 10 * log10(pwr);

% normalize psd_seg?
% psd_norm = psd ./ sum(psd,2);  % normalize by total power
% psd_seg = 10*log10(psd_seg);

% Visualize results
if vis
    figure('color', 'w');

    % % FFT power
    % subplot(4,1,1)
    % plot(f, pwr, 'k', 'LineWidth', 1);
    % box on; axis tight
    % ylabel('Power (µV^2/Hz)');
    % title('FFT');
    % 
    % % FFT power normalized to dB
    % subplot(4,1,2)
    % plot(f, pwr_norm, 'k', 'LineWidth', 1); hold on
    % ylabel('Power (dB)')
    % % xlabel('Frequency (Hz)')
    % box on; axis tight
    % title('FFT');

    % PSD
    % subplot(4,1,3)
    subplot(2,1,1)
    plot(f, psd, 'k', 'LineWidth', 1); hold on
    title('Power Spectral Density');
    ylabel('Power (µV^2/Hz)');
    box on; axis tight

    % PSD normalized to dB
    % subplot(4,1,4)
    subplot(2,1,2)
    plot(f, psd_norm, 'k', 'LineWidth', 1); hold on
    title('Power Spectral Density');
    ylabel('Power (dB)')
    xlabel('Frequency (Hz)')
    box on; axis tight

    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','normal');
end



