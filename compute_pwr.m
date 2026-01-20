function [pwr, psd, f, psd_seg, t_seg, coh_seg, pcoh_seg, gpdc_seg, coh_f, band_pwr] = compute_pwr(data, fs, varargin)
% COMPUTE_PWR  Welch or multitaper PSD with optional normalization, aperiodic removal,
%              coherence, and band power timecourses.
%
% OUTPUTS:
%   pwr        - Power per bin [nChan x nFreq] (μV²)
%   psd        - Power spectral density [nChan x nFreq] (μV²/Hz)
%   f          - Frequency vector (Hz)
%   psd_seg    - PSD per segment [nChan x nFreq x nSeg] (optional)
%   t_seg      - Time at end of each segment (s) (optional)
%   coh_seg    - Coherence per segment [nChan x nChan x nFreq x nSeg] (optional)
%   pcoh_seg   - Partial coherence per segment (optional)
%   gpdc_seg   - Generalized PDC per segment (optional)
%   coh_f      - Frequency vector for coherence (Hz) (optional)
%   band_pwr   - Struct with band power timecourses (optional)
%
% Name-Value options (case insensitive):
%   'Overlap'              - Fraction overlap in [0,1) (default 0.5)
%   'FRange'               - [fmin fmax] Hz (default [1 fs/2])
%   'WinLength'            - Window length in seconds (default 2)
%   'Method'               - 'welch' or 'multitaper' (default 'welch')
%   'NFFT'                 - FFT length (default winSize)
%   'NW'                   - Time-bandwidth product for multitaper (default 2.5)
%   'Tapers'               - Number of tapers for multitaper (default 2*NW-1)
%   'Combine'              - 'median' or 'mean' (default 'median')
%   'Progress'             - true/false/1/0 enable progressbar (default auto)
%   'dBNormalize'          - true/false/1/0 convert to dB (10*log10) (default false)
%   'AperiodicMethod'      - 'fooof', 'robust', 'ols', or 'none' (default 'none')
%   'AperiodicExcludeBands'- Nx2 array [fmin fmax] Hz to exclude from fit (default [])
%   'AperiodicFminFit'     - Lower bound for fit (default 1)
%   'AperiodicFmaxFit'     - Upper bound for fit (default max(FRange))
%   'FooofSettings'        - Struct with FOOOF settings (default [])
%   'CompCoherence'        - true/false/1/0 compute MVAR coherence (default false)
%   'CompBandPower'        - true/false/1/0 compute band power timecourses (default false)
%   'BandDefs'             - Struct with band definitions (default: delta/theta/alpha/beta/gamma)
%   'PlotPSD'              - true/false/1/0 plot PSD before and after processing
%   'UseParallel'          - true/false/1/0 parallelize (default false)

%% Parse inputs
p = inputParser;
p.KeepUnmatched = false;
p.CaseSensitive = false;

addParameter(p, 'Overlap', 0.5, @(x) x >= 0 && x < 1);
addParameter(p, 'FRange', []);
addParameter(p, 'WinLength', 2, @(x) x > 0);
addParameter(p, 'Method', 'welch', @(x) ismember(lower(x), {'welch', 'multitaper'}));
addParameter(p, 'NFFT', []);
addParameter(p, 'NW', 2.5, @(x) x > 0);
addParameter(p, 'Tapers', [], @(x) isempty(x) || (x > 0 && mod(x,1)==0));
addParameter(p, 'Combine', 'median', @(x) ismember(lower(x), {'median', 'mean'}));
addParameter(p, 'Progress', []);
addParameter(p, 'dBNormalize', false);
addParameter(p, 'AperiodicMethod', 'none', @(x) ismember(lower(x), {'fooof', 'robust', 'ols', 'none'}));
addParameter(p, 'AperiodicExcludeBands', []);
addParameter(p, 'AperiodicFminFit', []);
addParameter(p, 'AperiodicFmaxFit', []);
addParameter(p, 'FooofSettings', struct());
addParameter(p, 'CompCoherence', false);
addParameter(p, 'CompBandPower', false);
addParameter(p, 'BandDefs', []);
addParameter(p, 'PlotPSD', false);
addParameter(p, 'UseParallel', false);

parse(p, varargin{:});
opts = p.Results;

% Convert numeric 0/1 to logical for all boolean parameters
opts.Progress = to_logical(opts.Progress, []);
opts.dBNormalize = to_logical(opts.dBNormalize, false);
opts.CompCoherence = to_logical(opts.CompCoherence, false);
opts.CompBandPower = to_logical(opts.CompBandPower, false);
opts.PlotPSD = to_logical(opts.PlotPSD, false);
opts.UseParallel = to_logical(opts.UseParallel, false);

% Initialize optional outputs
coh_seg = []; pcoh_seg = []; gpdc_seg = []; coh_f = []; 
t_seg = []; psd_seg = []; band_pwr = struct();

% Progress settings (disabled for parallel mode)
if isempty(opts.Progress)
    opts.Progress = exist('progressbar', 'file') == 2;
end
usePar = opts.UseParallel && ~isempty(ver('parallel'));
if usePar && isempty(gcp('nocreate'))
    parpool;
end
usePar = usePar && ~isempty(gcp('nocreate'));
if usePar
    usePB = false;  % Disable progressbar in parallel mode
else
    usePB = opts.Progress;
end

%% Validate and reshape data
if size(data, 1) > size(data, 2)
    data = data';
end
if ndims(data) > 2
    data = reshape(data, size(data, 1), []);
end
[nChan, nSamples] = size(data);

if nChan > nSamples
    warning('compute_pwr:orientation', ...
        'Data appears to have more channels (%d) than samples (%d). Consider transposing.', ...
        nChan, nSamples);
end

% Frequency range
if isempty(opts.FRange)
    opts.FRange = [1 fs/2];
end
opts.FRange(1) = max(opts.FRange(1), 0);
opts.FRange(2) = min(opts.FRange(2), fs/2);

% Band definitions
if isempty(opts.BandDefs)
    opts.BandDefs = struct(...
        'delta', [0.5 4], ...
        'theta', [4 8], ...
        'alpha', [8 13], ...
        'beta', [13 30], ...
        'gamma', [30 80]);
end

%% Compute PSD using fixed windows
winSize = max(1, round(fs * opts.WinLength));
if winSize > nSamples
    winSize = nSamples;
    warning('compute_pwr:winSize', ...
        'Window size exceeds data length. Using full data length: %.2f s', winSize/fs);
end

step = max(1, winSize - floor(opts.Overlap * winSize));
segIdx = 1:step:(nSamples - winSize + 1);
if isempty(segIdx)
    segIdx = 1;
end
nSeg = numel(segIdx);

% FFT parameters
if isempty(opts.NFFT)
    nfft = winSize;
else
    nfft = opts.NFFT;
end
K = floor(nfft/2) + 1;
f = (0:(K-1)) * (fs/nfft);
if K > 1
    df = f(2) - f(1);
else
    df = fs/nfft;
end

% Prepare tapers
if strcmpi(opts.Method, 'multitaper')
    if isempty(opts.Tapers)
        Kt = max(1, floor(2*opts.NW - 1));
    else
        Kt = opts.Tapers;
    end
    [tapers, ~] = dpss(winSize, opts.NW, Kt);
else
    tapers = hamming(winSize);
end

%% Compute PSD per segment
psd_seg = nan(nChan, K, nSeg);
t_seg = nan(nSeg, 1);

if usePar
    fprintf('Computing PSD (parallel, %d segments)...\n', nSeg);

    dq = parallel.pool.DataQueue;
    nDone = 0;
    if opts.Progress
        afterEach(dq, @update_psd_progress);
    end

    parfor iSeg = 1:nSeg
        [psd_i, t_i] = psd_one_segment(data, segIdx(iSeg), winSize, fs, ...
            nfft, K, opts.Method, tapers);
        psd_seg(:, :, iSeg) = psd_i;
        t_seg(iSeg) = t_i;

        if opts.Progress
            send(dq, iSeg);
        end
    end

    if opts.Progress
        fprintf('PSD computation complete.\n');
    end
else
    progress_init('PSD computation', nSeg, usePB);
    for iSeg = 1:nSeg
        [psd_i, t_i] = psd_one_segment(data, segIdx(iSeg), winSize, fs, ...
            nfft, K, opts.Method, tapers);
        psd_seg(:, :, iSeg) = psd_i;
        t_seg(iSeg) = t_i;
        progress_step(iSeg, nSeg, usePB);
    end
    progress_finish(usePB);
end

%% Combine segments
psd = combine_segments(psd_seg, 3, opts.Combine);

%% Crop to frequency range
mask = f >= opts.FRange(1) & f <= opts.FRange(2);
f = f(mask);
psd = psd(:, mask);
psd_seg = psd_seg(:, mask, :);

% Update df for cropped range
if sum(mask) > 1
    df = mean(diff(f));
end

%% Store original PSD for plotting
psd_raw = psd;

%% Aperiodic correction
if ~strcmpi(opts.AperiodicMethod, 'none')
    % Set fit range
    if isempty(opts.AperiodicFminFit)
        fmin_fit = max(1, opts.FRange(1));
    else
        fmin_fit = max(opts.AperiodicFminFit, opts.FRange(1));
    end
    
    if isempty(opts.AperiodicFmaxFit)
        fmax_fit = opts.FRange(2);
    else
        fmax_fit = min(opts.AperiodicFmaxFit, opts.FRange(2));
    end
    
    % Create fit mask
    fitMask = (f >= max(fmin_fit, eps)) & (f <= fmax_fit);
    
    % Exclude specified bands
    if ~isempty(opts.AperiodicExcludeBands)
        for ii = 1:size(opts.AperiodicExcludeBands, 1)
            fitMask = fitMask & ~(f >= opts.AperiodicExcludeBands(ii,1) & ...
                f <= opts.AperiodicExcludeBands(ii,2));
        end
    end
    
    % Fit and remove aperiodic component
    if strcmpi(opts.AperiodicMethod, 'fooof')
        psd = remove_aperiodic_fooof(psd, f, fitMask, opts.FooofSettings);
    else
        psd = remove_aperiodic_regression(psd, f, fitMask, opts.AperiodicMethod);
    end
end

%% dB normalization
if opts.dBNormalize
    % Store linear for plotting
    psd_linear = psd;
    % Convert to dB (handle zeros/negatives)
    psd(psd <= 0) = eps;
    psd = 10 * log10(psd);
end

%% Compute power per bin
pwr = psd * df;

%% Coherence computation
if opts.CompCoherence
    if ~exist('fdMVAR_5order', 'file')
        warning('compute_pwr:noMVAR', ...
            'fdMVAR_5order not found. Skipping coherence computation.');
    else
        [coh_seg, pcoh_seg, gpdc_seg, coh_f] = compute_coherence(...
            data, fs, segIdx, winSize, nfft, opts.FRange, usePB, usePar, nChan);
    end
end

%% Band power timecourses
if opts.CompBandPower
    band_pwr = compute_band_power(psd_seg, f, t_seg, opts.BandDefs, opts.Combine);
end

%% Visualization
if opts.PlotPSD
    plot_psd(f, psd_raw, psd, opts.AperiodicMethod, opts.dBNormalize);
    
    if opts.CompBandPower && ~isempty(band_pwr)
        plot_band_timecourses(band_pwr);
    end
end

end

%% ======================== HELPER FUNCTIONS ========================

function [psd_i, t_end] = psd_one_segment(x, startIdx, winSize, fs, nfft, K, method, tapers)
% Compute PSD for one segment using proper Welch/multitaper scaling
idx = startIdx:(startIdx + winSize - 1);
t_end = idx(end) / fs;
xseg = x(:, idx);

if strcmpi(method, 'multitaper')
    % Multitaper method
    Sacc = zeros(size(xseg, 1), K);
    for kT = 1:size(tapers, 2)
        w = tapers(:, kT);
        % Compute scaling factor
        U = sum(w.^2);  % Energy of taper
        
        % Apply taper and compute FFT
        Xw = xseg .* w.';
        Y = fft(Xw, nfft, 2);
        Yp = Y(:, 1:K);
        
        % Compute PSD for this taper: S(f) = (1/(fs*U)) * |X(f)|^2
        Sxx = (1 / (fs * U)) * abs(Yp).^2;
        
        % Convert two-sided to single-sided (except DC and Nyquist)
        if rem(nfft, 2) == 0
            % Even nfft: double all except DC (1) and Nyquist (K)
            Sxx(:, 2:K-1) = 2 * Sxx(:, 2:K-1);
        else
            % Odd nfft: double all except DC (1)
            Sxx(:, 2:K) = 2 * Sxx(:, 2:K);
        end
        
        Sacc = Sacc + Sxx;
    end
    % Average across tapers
    psd_i = Sacc / size(tapers, 2);
else
    % Welch method with Hamming window
    w = tapers(:);
    % Compute scaling factor
    U = sum(w.^2);  % Energy of window
    
    % Apply window and compute FFT
    Xw = xseg .* w.';
    Y = fft(Xw, nfft, 2);
    Yp = Y(:, 1:K);
    
    % Compute PSD: S(f) = (1/(fs*U)) * |X(f)|^2
    psd_i = (1 / (fs * U)) * abs(Yp).^2;
    
    % Convert two-sided to single-sided
    if rem(nfft, 2) == 0
        psd_i(:, 2:K-1) = 2 * psd_i(:, 2:K-1);
    else
        psd_i(:, 2:K) = 2 * psd_i(:, 2:K);
    end
end
end

function A = combine_segments(X, dim, how)
% Combine segments using median or mean
switch lower(how)
    case 'median'
        A = median(X, dim, 'omitnan');
    case 'mean'
        A = mean(X, dim, 'omitnan');
    otherwise
        error('Combine must be ''median'' or ''mean''.');
end
end

function psd_corr = remove_aperiodic_fooof(psd, f, fitMask, fooof_settings)
% Remove aperiodic component using FOOOF
% This is a placeholder - actual FOOOF implementation would be needed

if ~exist('fooof', 'file')
    warning('compute_pwr:noFOOOF', ...
        'FOOOF not found. Falling back to OLS regression.');
    psd_corr = remove_aperiodic_regression(psd, f, fitMask, 'ols');
    return;
end

nChan = size(psd, 1);
psd_corr = psd;

% Default FOOOF settings
freq_range = [min(f(fitMask)) max(f(fitMask))];
if ~isfield(fooof_settings, 'peak_width_limits')
    fooof_settings.peak_width_limits = [0.5 12];
end
if ~isfield(fooof_settings, 'max_n_peaks')
    fooof_settings.max_n_peaks = Inf;
end
if ~isfield(fooof_settings, 'min_peak_height')
    fooof_settings.min_peak_height = 0;
end

for c = 1:nChan
    try
        % Run FOOOF
        fooof_results = fooof(f, psd(c, :), freq_range, fooof_settings);
        
        % Subtract aperiodic component
        ap_fit = fooof_results.ap_fit;
        psd_corr(c, :) = psd(c, :) - ap_fit;
        psd_corr(c, psd_corr(c, :) < 0) = 0;
    catch ME
        warning('compute_pwr:fooofFail', ...
            'FOOOF failed for channel %d: %s. Using OLS.', c, ME.message);
        psd_corr(c, :) = remove_aperiodic_regression(psd(c, :), f, fitMask, 'ols');
    end
end
end



function psd_corr = remove_aperiodic_regression(psd, f, fitMask, method)
% Remove aperiodic component using regression (OLS or robust)
nChan = size(psd, 1);
psd_fit = nan(size(psd));

logf_fit = log10(f(fitMask));
X = [ones(numel(logf_fit), 1) logf_fit(:)];

for c = 1:nChan
    y = psd(c, fitMask);
    y(y <= 0) = eps;
    logy = log10(y(:));
    
    if strcmpi(method, 'robust')
        % Iteratively reweighted least squares
        b = X \ logy;
        for iter = 1:5
            r = logy - X * b;
            s = 1.4826 * mad(r, 1);
            if s == 0
                s = eps;
            end
            w = huber_weights(r / s);
            W = diag(w);
            b = (X' * W * X) \ (X' * W * logy);
        end
    else
        % OLS
        b = X \ logy;
    end
    
    % Compute fit over entire frequency range
    logSfit = b(1) + b(2) * log10(f);
    psd_fit(c, :) = 10.^logSfit;
end

% Subtract aperiodic component
psd_corr = psd - psd_fit;
psd_corr(psd_corr < 0) = 0;
end

function w = huber_weights(r)
% Huber weighting function
c = 1.345;
a = abs(r);
w = ones(size(r));
w(a > c) = c ./ a(a > c);
end

function [coh_seg, pcoh_seg, gpdc_seg, coh_f] = compute_coherence(...
    data, fs, segIdx, winSize, nfft, fRange, usePB, usePar, nChan)
% Compute MVAR-based coherence measures per segment

nSeg = numel(segIdx);
Kc = floor(nfft/2) + 1;

coh_seg = nan(nChan, nChan, Kc, nSeg);
pcoh_seg = nan(nChan, nChan, Kc, nSeg);
gpdc_seg = nan(nChan, nChan, Kc, nSeg);
Su = eye(nChan);

% Get frequency axis
idx0 = segIdx(1):(segIdx(1) + winSize - 1);
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, cf0] = fdMVAR_5order(data(:, idx0), Su, nfft, fs);
coh_f = cf0(:).';

if usePar
    fprintf('Computing coherence (parallel, %d segments)...\n', nSeg);

    dq = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(dq, @update_coh_progress);

    parfor iSeg = 1:nSeg
        idx = segIdx(iSeg):(segIdx(iSeg) + winSize - 1);
        xw = data(:, idx);
        [~, ~, ~, gpdc, ~, coh, pcoh] = fdMVAR_5order(xw, Su, nfft, fs);
        coh_seg(:, :, :, iSeg) = abs(coh);
        pcoh_seg(:, :, :, iSeg) = abs(pcoh);
        gpdc_seg(:, :, :, iSeg) = abs(gpdc);

        send(dq, iSeg);
    end

    fprintf('Coherence computation complete.\n');
else
    progress_init('Coherence', nSeg, usePB);
    for iSeg = 1:nSeg
        idx = segIdx(iSeg):(segIdx(iSeg) + winSize - 1);
        xw = data(:, idx);
        [~, ~, ~, gpdc, ~, coh, pcoh] = fdMVAR_5order(xw, Su, nfft, fs);
        coh_seg(:, :, :, iSeg) = abs(coh);
        pcoh_seg(:, :, :, iSeg) = abs(pcoh);
        gpdc_seg(:, :, :, iSeg) = abs(gpdc);
        progress_step(iSeg, nSeg, usePB);
    end
    progress_finish(usePB);
end

% Crop to frequency range
if ~isempty(coh_f)
    cmask = coh_f >= fRange(1) & coh_f <= fRange(2);
    coh_f = coh_f(cmask);
    coh_seg = coh_seg(:, :, cmask, :);
    pcoh_seg = pcoh_seg(:, :, cmask, :);
    gpdc_seg = gpdc_seg(:, :, cmask, :);
end

function update_coh_progress(~)
    nDone = nDone + 1;
    if mod(nDone, max(1, round(nSeg/10))) == 0 || nDone == nSeg
        fprintf('  Coherence: %3.0f%% complete\n', 100*nDone/nSeg);
    end
end

end

function band_pwr = compute_band_power(psd_seg, f, t_seg, band_defs, combine_method)
% Compute power timecourses for each frequency band
bands = fieldnames(band_defs);
band_pwr = struct();

for b = 1:length(bands)
    band_name = bands{b};
    band_range = band_defs.(band_name);
    
    % Find frequencies in this band
    idx = f >= band_range(1) & f < band_range(2);
    
    if ~any(idx)
        warning('compute_pwr:emptyBand', ...
            'No frequencies found in %s band [%.1f-%.1f Hz]', ...
            band_name, band_range(1), band_range(2));
        continue;
    end
    
    % Average power across channels and frequencies
    band_data = psd_seg(:, idx, :);
    switch lower(combine_method)
        case 'median'
            timecourse = squeeze(median(median(band_data, 2, 'omitnan'), 1, 'omitnan'));
        case 'mean'
            timecourse = squeeze(mean(mean(band_data, 2, 'omitnan'), 1, 'omitnan'));
    end
    
    band_pwr.(band_name) = struct(...
        'timecourse', timecourse, ...
        't', t_seg, ...
        'freq_range', band_range);
end
end

function plot_psd(f, psd_raw, psd_processed, aperiodic_method, db_normalize)
% Plot PSD before and after processing
figure('Color', 'w', 'Position', [100 100 1000 500]);

subplot(1, 2, 1);
plot(f, psd_raw.', 'LineWidth', 1.5);
title('Original PSD', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 11);
if db_normalize
    ylabel('Power (dB)', 'FontSize', 11);
else
    ylabel('Power (\muV^2/Hz)', 'FontSize', 11);
end
box on;
set(gca, 'FontSize', 10);

subplot(1, 2, 2);
plot(f, psd_processed.', 'LineWidth', 1.5);
if ~strcmpi(aperiodic_method, 'none') && db_normalize
    title_str = 'Aperiodic-Corrected & dB-Normalized PSD';
elseif ~strcmpi(aperiodic_method, 'none')
    title_str = 'Aperiodic-Corrected PSD';
elseif db_normalize
    title_str = 'dB-Normalized PSD';
else
    title_str = 'Processed PSD';
end
title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 11);
if db_normalize
    ylabel('Power (dB)', 'FontSize', 11);
else
    ylabel('Power (\muV^2/Hz)', 'FontSize', 11);
end
box on;
set(gca, 'FontSize', 10);
end

function plot_band_timecourses(band_pwr)
% Plot frequency band power over time
bands = fieldnames(band_pwr);
nBands = length(bands);

figure('Color', 'w', 'Position', [100 100 1000 800]);

for b = 1:nBands
    subplot(nBands, 1, b);
    hold on;
    
    data = band_pwr.(bands{b});
    t = data.t;
    tc = data.timecourse;
    
    % Plot raw and smoothed
    plot(t, tc, '.', 'MarkerSize', 8, 'Color', [0.7 0.7 0.7]);
    plot(t, smooth(tc, min(30, length(tc))), 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    
    title(sprintf('%s (%.1f-%.1f Hz)', upper(bands{b}), ...
        data.freq_range(1), data.freq_range(2)), ...
        'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Power (\muV^2/Hz)', 'FontSize', 10);
    
    if b == nBands
        xlabel('Time (s)', 'FontSize', 10);
    end
    
    box on;
    set(gca, 'FontSize', 9);
end
end

function val = to_logical(input, default_val)
% Convert 0/1 or logical to logical, with default
if isempty(input)
    if isempty(default_val)
        val = [];
    else
        val = logical(default_val);
    end
elseif isnumeric(input)
    val = logical(input);
elseif islogical(input)
    val = input;
else
    error('Parameter must be logical or numeric 0/1');
end
end

function tf = hasParallel()
% Check if Parallel Computing Toolbox is available and pool exists
tf = ~isempty(ver('parallel'));
if tf
    pool = gcp('nocreate');
    tf = ~isempty(pool);
    if ~tf
        warning('compute_pwr:noPool', ...
            'Parallel Computing Toolbox available but no pool started. Use parpool() or set UseParallel=false.');
    end
end
end

%% Progress tracking functions
function progress_init(label, nTotal, usePB)
persistent p_total p_count p_label p_usePB
p_total = max(1, nTotal);
p_count = 0;
p_label = label;
p_usePB = usePB;

if p_usePB
    progressbar(0);
end
end

function progress_step(i, n, usePB)
if usePB
    progressbar(i / n);
end
end

function progress_finish(usePB)
if usePB
    progressbar(1);
end
end


function update_psd_progress(~)
    nDone = nDone + 1;
    if mod(nDone, max(1, round(nSeg/10))) == 0 || nDone == nSeg
        fprintf('  PSD: %3.0f%% complete\n', 100*nDone/nSeg);
    end
end
