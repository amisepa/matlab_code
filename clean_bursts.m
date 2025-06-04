function bad_mask = clean_bursts(signal, fs, window_len, stepsize, zthreshold)
% CLEAN_BURSTS: Detects RMS bursts and flatline periods in real-time.
%
% Inputs:
%   signal      - [channels x samples] matrix
%   fs          - Sampling rate in Hz
%   window_len  - Window length in seconds
%   stepsize    - Step size in samples
%   zthreshold  - Z-score threshold for RMS bursts
%
% Output:
%   bad_mask    - Logical [channels x samples] matrix of flagged samples

if nargin < 5 || isempty(zthreshold)
    zthreshold = 4;
end

[chan, num_samples] = size(signal);
samples_per_window = round(fs * window_len);
offsets = 1:stepsize:(num_samples - samples_per_window + 1);

% Precompute window indices
win_idx = bsxfun(@plus, (0:samples_per_window - 1)', offsets);  % [samples_per_window x num_windows]

% Initialize mask
bad_mask = false(chan, num_samples);

% Channel-wise burst and flatline detection
for ch = 1:chan
    x = signal(ch, :);
    xw = x(win_idx);  % windowed signal [samples x windows]

    % RMS burst detection
    rms_vals = sqrt(mean(xw.^2, 1));  % 1 x num_windows
    mu = median(rms_vals);
    sigma = 1.4826 * mad(rms_vals, 1);
    z = abs((rms_vals - mu) / sigma);
    burst_wins = find(z > zthreshold);

    % Flatline detection (based on std per window)
    std_vals = std(xw, 0, 1);
    flat_thresh = median(std_vals) * 0.2;
    % flat_thresh = 1e-3;  % for signals in ~µV or mV range
    flat_wins = find(std_vals < flat_thresh);

    % Combine window indices
    bad_wins = union(burst_wins, flat_wins);
    if ~isempty(bad_wins)
        idx = win_idx(:, bad_wins);
        idx = idx(:);
        idx(idx > num_samples) = [];
        bad_mask(ch, idx) = true;
    end
end
end
