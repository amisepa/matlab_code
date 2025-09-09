function [metrics, win_time] = compute_eeg_signal_features(signal, fs, win_len, overlap)
% compute_eeg_signal_features
% Sliding window EEG quality and morphology features.
%
% Inputs
%   signal   1D EEG vector
%   fs       sampling rate in Hz
%   win_len  window length in seconds  default 2
%   overlap  logical for 50 percent overlap  default false
%
% Outputs
%   metrics  struct of per window features  fields are vectors
%   win_time window start times in minutes
%
% Example
%   fs = 250;
%   t = 0:1/fs:60;
%   x = 20*sin(2*pi*10*t) + 0.5*randn(size(t));
%   [M, wt] = compute_eeg_signal_features(x, fs, 2, true);

if nargin < 3 || isempty(win_len), win_len = 2; end
if nargin < 4 || isempty(overlap), overlap = false; end

signal = signal(:);                         % force column
if size(signal,2) ~= 1
    error('Input signal must be 1D.');
end

N = numel(signal);
wlen = max(1, round(win_len * fs));
step = overlap * floor(wlen/2) + (~overlap) * wlen;
seg_idx = 0:step:(N - wlen);
seg_idx = seg_idx(:);

% Predefine feature containers
keys = { ...
    'flat', 'pSNR', 'tSNR', 'SNR1', 'SNR2', 'residual_mad', ...
    'rms', 'kurt', 'skew', 'shannon_entropy', 'spec_entropy', ...
    'mean_slope', 'slope_var', 'zc', 'zc_int', 'fractal', ...
    'lag1_autocorr', 'signal_energy', 'signal_power' ...
    };

for k = 1:numel(keys)
    metrics.(keys{k}) = [];
end

% EEG filters
% High pass 0.75 Hz, Band pass 0.75 to 45 Hz
try
    [bhp, ahp] = butter(2, 0.75/(fs/2), 'high');
    [bbp, abp] = butter(2, [0.75 45]/(fs/2), 'bandpass');
    hp_filt = filtfilt(bhp, ahp, double(signal));
    bp_filt = filtfilt(bbp, abp, double(signal));
catch ME
    error('Filter failure: %s', ME.message);
end

% FFT frequency vector for each window length
freqs = (0:(wlen/2 - 1))' * (fs / wlen);

% EEG signal band mask 4 to 21 Hz
eeg_mask = (freqs >= 4) & (freqs <= 21);

% Hann window for spectral calc
win = hann(wlen, 'periodic');

for ii = 1:numel(seg_idx)
    startIdx = seg_idx(ii) + 1;
    stopIdx  = startIdx + wlen - 1;
    if stopIdx > N
        break
    end

    hp_seg = hp_filt(startIdx:stopIdx);
    bp_seg = bp_filt(startIdx:stopIdx);

    if numel(hp_seg) < 8
        continue
    end

    % PSD via windowed FFT
    xw = hp_seg .* win;
    X  = fft(xw);
    P2 = abs(X).^2 / (fs * wlen);
    P1 = P2(1:floor(wlen/2));
    P1 = P1 + eps;
    psd_norm = P1 / sum(P1);
    
    sig_pow   = sum(P1(eeg_mask));
    noise_pow = sum(P1(~eeg_mask) + eps);

    % flat percent
    metrics.flat(end+1,1) = sum(abs(diff(hp_seg)) < 1e-6) / wlen * 100;

    % power SNR
    metrics.pSNR(end+1,1) = 10 * safe_log10(sig_pow / noise_pow);

    % time SNR  AC amplitude over HF noise after trend removal
    ac_amp = max(bp_seg) - min(bp_seg);
    kernel_len = min(max(3, round(0.2 * fs)), floor(numel(bp_seg)/2));
    trend = movmean(bp_seg, kernel_len);
    rms_noise = sqrt(mean((bp_seg - trend).^2));
    metrics.tSNR(end+1,1) = 20 * safe_log10(ac_amp / (rms_noise + eps));

    % variance SNR forms
    metrics.SNR1(end+1,1) = 10 * safe_log10(var(bp_seg) / (var(hp_seg - bp_seg) + eps));
    metrics.SNR2(end+1,1) = 10 * safe_log10(sig_pow / noise_pow);

    % residual MAD between high pass and band pass
    metrics.residual_mad(end+1,1) = median(abs(hp_seg - bp_seg));

    % amplitude stats
    metrics.rms(end+1,1)  = sqrt(mean(bp_seg.^2));
    metrics.kurt(end+1,1) = kurtosis(bp_seg);
    metrics.skew(end+1,1) = skewness(bp_seg);

    % Shannon entropy from amplitude histogram
    [hcounts, ~] = histcounts(bp_seg, 10, 'Normalization', 'probability');
    hcounts = hcounts(hcounts > 0);
    metrics.shannon_entropy(end+1,1) = -sum(hcounts .* log2(hcounts));

    % Spectral entropy
    metrics.spec_entropy(end+1,1) = -sum(psd_norm .* log2(psd_norm));

    % slope stats
    slopes = diff(bp_seg);
    metrics.mean_slope(end+1,1) = mean(slopes);
    metrics.slope_var(end+1,1)  = var(slopes);

    % zero crossings and intervals
    zc_idx = find(diff(sign(bp_seg)) ~= 0);
    metrics.zc(end+1,1) = numel(zc_idx);
    metrics.zc_int(end+1,1) = ternary(numel(zc_idx) > 1, mean(diff(zc_idx))/fs, NaN);

    % Petrosian fractal dimension approximation
    n = numel(bp_seg);
    Nz = max(numel(zc_idx), 1);
    metrics.fractal(end+1,1) = log10(n) / (log10(n) + log10(n/(n + 0.4*Nz)));

    % lag 1 autocorrelation
    acf = xcorr(bp_seg, 1, 'unbiased');
    % xcorr returns [lag -1, lag 0, lag +1] for maxlag 1
    r0 = acf(2);
    r1 = acf(3);
    if r0 ~= 0
        metrics.lag1_autocorr(end+1,1) = r1 / r0;
    else
        metrics.lag1_autocorr(end+1,1) = NaN;
    end

    % energy and power
    metrics.signal_energy(end+1,1) = sum(bp_seg.^2);
    metrics.signal_power(end+1,1)  = mean(bp_seg.^2);
end

win_time = seg_idx / fs / 60;  % minutes
end

function y = safe_log10(x)
y = log10(max(x, realmin));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
