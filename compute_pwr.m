function [pwr, pwr_norm, psd, psd_norm, f, psd_seg, t_seg, coh_seg, pcoh_seg, gpdc_seg, coh_f] = compute_pwr(data, fs, opts)
% COMPUTE_PWR  Welch or multitaper PSD with correct scaling, robust combine, optional coherence,
%               optional multiresolution windows, and optional 1/f aperiodic fit removal.
%
% Copyright (C) Cedric Cannard 2021-2025
%
% Usage
%   [pwr, pwr_norm, psd, psd_norm, f, psd_seg, t_seg, coh_seg, pcoh_seg, gpdc_seg, coh_f] = ...
%       compute_pwr(data, fs, opts)
%
% Inputs
%   data    [nChan x nSamples] or [nChan x nSamples x nTrials]
%   fs      sampling rate in Hz
%   opts    struct of options fields below. Any field can be omitted.
%
% Options in opts
%   overlap           fraction overlap in [0,1), default 0.5
%   fRange            [fmin fmax] Hz, default [1 fs/2]
%   winLength         window length in seconds, default 2
%   method            'welch' or 'multitaper', default 'welch'
%   nfft              FFT length, default winSize
%   detrend           'mean' (default) | 'linear' | 'none'
%   trimPercent       20 - percent trimmed from each tail for trimmean
%   combine           'trimmean' (default) | 'median' | 'mean'
%   progress          true to show progressbar if available, default auto-detect
%   timeBW            DPSS time-halfbandwidth NW for multitaper, default 2.5
%   numTapers         number of DPSS tapers K, default floor(2*NW - 1)
%   adaptive_bands    struct array with fields .fmin .fmax .cycles - stitches multiwindow PSD
%   comp_coherence    true to compute MVAR coherence per segment, default false
%   aperiodic.enable        true to compute 1/f removal, default true
%   aperiodic.method        'ols' (default) or 'robust'
%   aperiodic.excludeBands  Nx2 array of [fmin fmax] Hz to exclude from fit, default []
%   aperiodic.fmin_fit      lower bound for fit, default max(1, fRange(1))
%   aperiodic.fmax_fit      upper bound for fit, default fRange(2)
%   plot_psd          true to plot PSD before and after 1/f removal
%   plot_bands        true to plot delta/theta/alpha/beta timecourses from psd_seg
%
% Outputs
%   psd, pwr          raw one-sided spectra in linear units
%   psd_norm, pwr_norm  oscillation-only (aperiodic-removed) in linear units
%   psd_seg           segment-wise PSD [nChan x nFreq x nSegments] from the global window
%   t_seg             segment end times in seconds
%   coh_seg, pcoh_seg, gpdc_seg, coh_f  coherence products if comp_coherence = true
%
% Examples
%   % Minimal, Welch, 2 s windows, default overlap, aperiodic removal on
%   opts = struct('plot_psd',true);
%   [pwr, pwr_osc, psd, psd_osc, f] = compute_pwr(data, fs, opts);
%
%   % Multitaper with NW=3, K=5, linear detrend, wider range, band plots
%   opts = struct('method','multitaper','timeBW',3,'numTapers',5, ...
%                 'detrend','linear','fRange',[1 45], ...
%                 'plot_psd',true,'plot_bands',true);
%   [pwr, pwr_osc, psd, psd_osc, f, psd_seg, t_seg] = compute_pwr(data, fs, opts);
%
%   % Adaptive windows stitched for classic EEG bands
%   % Adaptive windows mode computes separate PSDs with window length tailored
%   % to each frequency band (e.g., cycles per center frequency) and stitches
%   % those bins into a single spectrum to better capture band-limited oscillations.
%   bands = [struct('fmin',1,'fmax',3,'cycles',8), ...
%            struct('fmin',4,'fmax',7.5,'cycles',8), ...
%            struct('fmin',8,'fmax',13,'cycles',8), ...
%            struct('fmin',13,'fmax',30,'cycles',10)];
%   opts = struct('adaptive_bands',bands,'plot_psd',true);
%   [pwr, pwr_osc, psd, psd_osc, f] = compute_pwr(data, fs, opts);
%
%   % Adaptive windows changes the window length with frequency band to tune
%   % time-frequency tradeoffs, while multitaper uses multiple DPSS tapers within
%   % a fixed window to reduce variance and leakage; you can even combine them
%   % by running multitaper inside each band-specific window.
%   bands = [struct('fmin',4,'fmax',7.5,'cycles',8), struct('fmin',8,'fmax',13,'cycles',8)];
%   opts = struct('method','multitaper','timeBW',3,'numTapers',5,'adaptive_bands',bands,'plot_psd',true);
%   [pwr, pwr_osc, psd, psd_osc, f] = compute_pwr(data, fs, opts);

% -------- defaults and arg parsing --------
if nargin < 3 || isempty(opts), opts = struct(); end

def  = struct('overlap',0.5,'fRange',[],'winLength',2,'method','welch','nfft',[], ...
              'detrend','mean','trimPercent',20,'combine','trimmean','progress',[], ...
              'timeBW',2.5,'numTapers',[],'adaptive_bands',[], ...
              'comp_coherence',false, ...
              'aperiodic',struct(), 'plot_psd',false, 'plot_bands',false);
apdef = struct('enable',true,'method','ols','excludeBands',[], 'fmin_fit',[], 'fmax_fit',[]);
opts  = merge_structs(def, opts);
opts.aperiodic = merge_structs(apdef, opts.aperiodic);

if isempty(opts.fRange), opts.fRange = [1 fs/2]; end
if isempty(opts.progress), opts.progress = exist('progressbar','file') == 2; end

% -------- shape checks --------
coh_seg = []; pcoh_seg = []; gpdc_seg = []; coh_f = []; t_seg = []; psd_seg = [];
if size(data,1) > size(data,2), data = data'; end
if ndims(data) > 2, data = reshape(data, size(data,1), []); end
[nChan, nSamples] = size(data);

% -------- global-window estimate --------
[Sxx_seg_global, Pbin_seg_global, f_global, t_seg, df] = ...
    estimate_psd_fixedwin(data, fs, opts.winLength, opts.method, opts.nfft, opts.overlap, nSamples, nChan, opts);

% Use global segments for time plots; compute averaged spectra from global
psd_seg = Sxx_seg_global;                 % keep for time-course plots
f       = f_global;
psd     = combine_segments(psd_seg, 3, opts);      % averaged PSD from global segments
pwr     = combine_segments(Pbin_seg_global, 3, opts);

% -------- adaptive bands: stitch at averaged level only --------
useAdaptive = isstruct(opts.adaptive_bands) && ~isempty(opts.adaptive_bands);
if useAdaptive
    bands = opts.adaptive_bands;
    for b = 1:numel(bands)
        fmin = bands(b).fmin; fmax = bands(b).fmax; cyc = bands(b).cycles;
    if ~(fmin < fmax && cyc > 0), warning('Skipping band %d: invalid fields.', b); continue; end
        f0 = sqrt(fmin*fmax);
        winLen_b = cyc / max(f0, eps);

        % Compute spectra with the band-specific window
        [Sxx_b, Pbin_b, f_b, ~, ~] = estimate_psd_fixedwin(data, fs, winLen_b, opts.method, opts.nfft, opts.overlap, nSamples, nChan, opts);
        psd_b = combine_segments(Sxx_b,   3, opts);
        pwr_b = combine_segments(Pbin_b,  3, opts);

        % Map only bins inside this band from f_b onto global f and override
        mask_b = f_b >= fmin & f_b <= fmax;
        if ~any(mask_b), continue; end
        fb = f_b(mask_b);
        [~, idxGlobal] = arrayfun(@(fv)findClosest(f, fv), fb);
        psd(:, idxGlobal) = psd_b(:, mask_b);
        pwr(:, idxGlobal) = pwr_b(:, mask_b);
    end
end

% -------- crop frequency range --------
mask = f >= opts.fRange(1) & f <= opts.fRange(2);
f         = f(mask);
psd       = psd(:, mask);
pwr       = pwr(:, mask);
psd_seg   = psd_seg(:, mask, :);

% -------- aperiodic 1/f fit and removal in linear units --------
psd_norm = []; pwr_norm = []; psd_fit = [];
if opts.aperiodic.enable
    fmin_fit = max(1, opts.fRange(1)); if ~isempty(opts.aperiodic.fmin_fit), fmin_fit = max(opts.aperiodic.fmin_fit, opts.fRange(1)); end
    fmax_fit = opts.fRange(2);          if ~isempty(opts.aperiodic.fmax_fit), fmax_fit = min(opts.aperiodic.fmax_fit, opts.fRange(2)); end
    fitMask  = (f >= max(fmin_fit, eps)) & (f <= fmax_fit);
    if ~isempty(opts.aperiodic.excludeBands)
        ex = false(size(f));
        for ii = 1:size(opts.aperiodic.excludeBands,1)
            ex = ex | (f >= opts.aperiodic.excludeBands(ii,1) & f <= opts.aperiodic.excludeBands(ii,2));
        end
        fitMask = fitMask & ~ex;
    end
    fitMask = fitMask(:).';
    logf_fit = log10(f(fitMask));
    psd_fit = nan(size(psd), 'like', psd);

    for c = 1:nChan
        y = psd(c, fitMask); y(y<=0) = eps;
        logy = log10(y);
        X = [ones(numel(logf_fit),1) logf_fit(:)];
        if strcmpi(opts.aperiodic.method,'robust')
            b = X \ logy(:);
            for iter = 1:5
                r = logy(:) - X*b;
                s = 1.4826*mad(r,1); if s==0, s = eps; end
                w = huberWeights(r/s); W = diag(w);
                b = (X'*W*X) \ (X'*W*logy(:));
            end
        else
            b = X \ logy(:);
        end
        logSfit = b(1) + b(2)*log10(f);
        psd_fit(c,:) = 10.^logSfit;
    end

    psd_osc = psd - psd_fit; 
    psd_osc(psd_osc < 0) = 0;
    psd_norm = psd_osc;

    if numel(f) > 1
        df_local = mean(diff(f));
    else
        df_local = df;   % from estimate_psd_fixedwin
    end
    pwr_norm = psd_norm .* df_local;
end

% -------- coherence using the global window only --------
if opts.comp_coherence
    winSize = max(1, round(fs*opts.winLength));
    step    = max(1, winSize - floor(opts.overlap*winSize));
    segIdx  = 1:step:(nSamples - winSize + 1); if isempty(segIdx), segIdx = 1; end
    nSeg    = numel(segIdx);
    nfftC   = isempty(opts.nfft) * winSize + ~isempty(opts.nfft) * opts.nfft;
    Kc      = floor(nfftC/2) + 1;
    coh_seg  = nan(nChan, nChan, Kc, nSeg);
    pcoh_seg = nan(nChan, nChan, Kc, nSeg);
    gpdc_seg = nan(nChan, nChan, Kc, nSeg);
    % if opts.progress, progressbar('Computing coherence'); end
    if opts.progress, progressbar(0, nSeg); end
    Su = eye(nChan);
    for iSeg = 1:nSeg
        idx = segIdx(iSeg):(segIdx(iSeg)+winSize-1);
        xw  = data(:, idx);
        switch lower(opts.detrend)
            case 'linear'
                t = (0:size(xw,2)-1); t = t - mean(t); denom = sum(t.^2);
                for cc = 1:nChan
                    slope = (t * double(xw(cc,:))')/denom;
                    xw(cc,:) = xw(cc,:) - (slope*t);
                    xw(cc,:) = xw(cc,:) - mean(xw(cc,:));
                end
            case 'none'
            otherwise
                xw = xw - mean(xw,2);
        end
        [~, ~, ~, gpdc, ~, coh, pcoh, ~, ~, ~, ~, cf] = fdMVAR_5order(xw, Su, nfftC, fs);
        coh_seg(:,:,:,iSeg)  = abs(coh);
        pcoh_seg(:,:,:,iSeg) = abs(pcoh);
        gpdc_seg(:,:,:,iSeg) = abs(gpdc);
        if iSeg == 1, coh_f = cf(:).'; end
        % if opts.progress, progressbar(iSeg/nSeg); end
        if opts.progress, progressbar(iSeg, nSeg); end
    end
    if ~isempty(coh_f)
        cmask   = coh_f >= opts.fRange(1) & coh_f <= opts.fRange(2);
        coh_f   = coh_f(cmask);
        coh_seg = coh_seg(:,:,cmask,:);
        pcoh_seg= pcoh_seg(:,:,cmask,:);
        gpdc_seg= gpdc_seg(:,:,cmask,:);
    end
end

% -------- visualization via opts.plot_* --------
if isfield(opts,'plot_psd') && opts.plot_psd
    figure('color','w');
    subplot(2,1,1); hold on
    plot(f, psd.', 'LineWidth', 1);
    if ~isempty(psd_fit)
        plot(f, mean(psd_fit,1,'omitnan'), '--', 'LineWidth', 1);
    end
    title('Power Spectral Density - raw');
    ylabel('Power (\muV^2/Hz)'); box on; axis tight

    subplot(2,1,2); hold on
    if ~isempty(psd_norm)
        plot(f, psd_norm.', 'LineWidth', 1);
        title('Power Spectral Density - oscillation only');
    else
        plot(f, psd.', 'LineWidth', 1);
        title('Power Spectral Density');
    end
    ylabel('Power (\muV^2/Hz)'); xlabel('Frequency (Hz)'); box on; axis tight
    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','normal');
end

if isfield(opts,'plot_bands') && opts.plot_bands && ~isempty(psd_seg)
    figure('color','w')
    subplot(4,1,1); hold on
    idx = f <= 3;
    data_ts = squeeze(trimmean(trimmean(psd_seg(:,idx,:),20,2),20,1));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Delta'); ylabel('\muV^2/Hz'); box on

    subplot(4,1,2); hold on
    idx = f > 3 & f < 8;
    data_ts = squeeze(trimmean(trimmean(psd_seg(:,idx,:),20,2),20,1));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Theta'); ylabel('\muV^2/Hz'); box on

    subplot(4,1,3); hold on
    idx = f >= 8 & f < 13;
    data_ts = squeeze(trimmean(trimmean(psd_seg(:,idx,:),20,2),20,1));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Alpha'); ylabel('\muV^2/Hz'); box on

    subplot(4,1,4); hold on
    idx = f >= 13 & f < 30;
    data_ts = squeeze(trimmean(trimmean(psd_seg(:,idx,:),20,2),20,1));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Beta'); ylabel('\muV^2/Hz'); xlabel('Time (s)'); box on

    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','normal');
end
end


%% Helper functions - keep at bottom

function A = combine_segments(X, dim, opts)
switch lower(opts.combine)
    case 'median', A = median(X, dim, 'omitnan');
    case 'mean',   A = mean(X, dim, 'omitnan');
    otherwise,     A = trimmean(X, opts.trimPercent, dim);
end
end

function [Sxx_seg, Pbin_seg, fvec, t_ends, df] = estimate_psd_fixedwin( ...
        x, fsL, winSec, method, nfftUser, overlap, nSamples, nChan, opts)
% Standalone helper. All needed parameters are passed explicitly.

winSize = max(1, round(fsL*winSec));
if winSize > nSamples, winSize = nSamples; end

step    = max(1, winSize - floor(overlap*winSize));
segIdx  = 1:step:(nSamples - winSize + 1);
if isempty(segIdx), segIdx = 1; end
nSeg    = numel(segIdx);

% FFT params
if isempty(nfftUser), nfft = winSize; else, nfft = nfftUser; end
K    = floor(nfft/2) + 1;
fvec = (0:(K-1))*(fsL/nfft);
if K > 1
    df = fvec(2) - fvec(1);
else
    df = fsL/nfft;
end

% Preallocate
Sxx_seg = nan(nChan, K, nSeg, 'like', x);
Pbin_seg= nan(nChan, K, nSeg, 'like', x);
t_ends  = nan(nSeg,1);

% Tapers
if strcmpi(method,'multitaper')
    NW = opts.timeBW;
    if isempty(opts.numTapers)
        Kt = max(1, floor(2*NW - 1));
    else
        Kt = opts.numTapers;
    end
    [tapers, ~] = dpss(winSize, NW, Kt);   % [winSize x Kt]
else
    tapers = hamming(winSize);             % [winSize x 1]
end

% if opts.progress, progressbar(sprintf('PSD %s', method)); end
if opts.progress, progressbar(0, nSeg); end
for iSeg = 1:nSeg
    idx          = segIdx(iSeg):(segIdx(iSeg)+winSize-1);
    t_ends(iSeg) = idx(end)/fsL;

    xseg = x(:, idx);                      % [nChan x winSize]

    % Detrend
    switch lower(opts.detrend)
        case 'linear'
            t = (0:winSize-1); t = t - mean(t);
            denom = sum(t.^2);
            for cc = 1:nChan
                slope = (t * double(xseg(cc,:))')/denom;
                xseg(cc,:) = xseg(cc,:) - (slope*t);
                xseg(cc,:) = xseg(cc,:) - mean(xseg(cc,:));
            end
        case 'none'
            % do nothing
        otherwise
            xseg = xseg - mean(xseg,2);
    end

    if strcmpi(method,'multitaper')
        % Average across DPSS tapers with window-power correction
        Sacc = zeros(nChan, K, 'like', xseg);
        for kT = 1:size(tapers,2)
            w   = tapers(:,kT).';
            W2  = sum(w.^2);
            Xw  = xseg .* w;
            Y   = fft(Xw, nfft, 2);
            Yp  = Y(:, 1:K);
            Sxx = (1/(fsL*W2)) * abs(Yp).^2;
            if rem(nfft,2) == 0
                pos = 2:K-1;          % even nfft: leave DC and Nyquist
            else
                pos = 2:K;            % odd nfft: leave DC
            end
            Sxx(:,pos) = 2*Sxx(:,pos);
            Sacc = Sacc + Sxx;
        end
        Sxx = Sacc / size(tapers,2);
    else
        % Welch with Hamming
        w   = tapers(:).';
        W2  = sum(w.^2);
        Xw  = xseg .* w;
        Y   = fft(Xw, nfft, 2);
        Yp  = Y(:, 1:K);
        Sxx = (1/(fsL*W2)) * abs(Yp).^2;   % window-power corrected
        if rem(nfft,2) == 0
            pos = 2:K-1;                  % even nfft: leave DC and Nyquist
        else
            pos = 2:K;                    % odd nfft: leave DC
        end
        Sxx(:,pos) = 2*Sxx(:,pos);        % one-sided
    end

    Sxx_seg(:,:,iSeg) = Sxx;              % µV^2/Hz if input is µV
    Pbin_seg(:,:,iSeg)= Sxx * df;         % per-bin power, µV^2

    % if opts.progress, progressbar(iSeg/nSeg); end
    if opts.progress, progressbar(iSeg, nSeg); end
end
end

function [d, idx] = findClosest(vec, val)
[~, idx] = min(abs(vec - val)); d = vec(idx);
end

function w = huberWeights(r)
c = 1.345; a = abs(r); w = ones(size(r)); w(a > c) = c ./ a(a > c);
end

function S = merge_structs(base, add)
if isempty(add)
    S = base; 
    return
end
S = base;
f = fieldnames(add);
for ii = 1:numel(f)
    S.(f{ii}) = add.(f{ii}); 
end
end
