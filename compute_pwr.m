function [pwr, pwr_norm, psd, psd_norm, f, psd_seg, t_seg, coh_seg, pcoh_seg, gpdc_seg, coh_f] = compute_pwr(data, fs, varargin)
% COMPUTE_PWR  Welch or multitaper PSD with correct scaling, robust combine, optional coherence,
%               optional 1/f removal, optional parallel, and clean progress tracking.
%
% Name–Value options (case insensitive)
%   'Overlap'              fraction overlap in [0,1) (default 0.5)
%   'FRange'               [fmin fmax] Hz (default [1 fs/2])
%   'WinLength'            window length in seconds (default 2)
%   'Method'               'welch' or 'multitaper' (default 'welch')
%   'NFFT'                 FFT length (default winSize)
%   'Combine'              'median' or 'mean' (default 'median')
%   'Progress'             true to enable progressbar if available (default auto)
%   'CompCoherence'        true to compute MVAR coherence per segment (default false)
%   'AperiodicEnable'      true to compute 1/f removal (default true)
%   'AperiodicMethod'      'ols' or 'robust' (default 'ols')
%   'AperiodicExcludeBands'  Nx2 array [fmin fmax] Hz to exclude from fit (default [])
%   'AperiodicFminFit'     lower bound for fit (default max(1, FRange(1)))
%   'AperiodicFmaxFit'     upper bound for fit (default FRange(2))
%   'PlotPSD'              true to plot PSD before and after 1/f removal
%   'PlotBands'            true to plot delta/theta/alpha/beta timecourses from psd_seg
%   'UseParallel'          true to parallelize across segments for PSD and coherence (default false)

% -------- parse name–value args --------
p = inputParser; p.KeepUnmatched = false; p.CaseSensitive = false;
addParameter(p,'Overlap',0.5);
addParameter(p,'FRange',[]);
addParameter(p,'WinLength',2);
addParameter(p,'Method','welch');
addParameter(p,'NFFT',[]);
addParameter(p,'Combine','median');           % only 'median' or 'mean'
addParameter(p,'Progress',[]);
addParameter(p,'CompCoherence',false);
addParameter(p,'AperiodicEnable',true);
addParameter(p,'AperiodicMethod','ols');
addParameter(p,'AperiodicExcludeBands',[]);
addParameter(p,'AperiodicFminFit',[]);
addParameter(p,'AperiodicFmaxFit',[]);
addParameter(p,'PlotPSD',false);
addParameter(p,'PlotBands',false);
addParameter(p,'UseParallel',false);
parse(p,varargin{:});
opts = p.Results;

% progress default
if isempty(opts.Progress), opts.Progress = exist('progressbar','file') == 2; end
usePB  = logical(opts.Progress);
usePar = logical(opts.UseParallel) && hasParallel();

% -------- shape checks --------
coh_seg = []; pcoh_seg = []; gpdc_seg = []; coh_f = []; t_seg = []; psd_seg = [];
if size(data,1) > size(data,2), data = data'; end
if ndims(data) > 2, data = reshape(data, size(data,1), []); end
[nChan, nSamples] = size(data);

if isempty(opts.FRange), opts.FRange = [1 fs/2]; end

% -------- global-window estimate --------
[Sxx_seg_global, Pbin_seg_global, f_global, t_seg, df] = ...
    estimate_psd_fixedwin(data, fs, opts.WinLength, opts.Method, opts.NFFT, opts.Overlap, nSamples, nChan, usePB, usePar, 'PSD');

psd_seg = Sxx_seg_global;
f       = f_global;
psd     = combine_segments(psd_seg, 3, opts.Combine);
pwr     = combine_segments(Pbin_seg_global, 3, opts.Combine);

% -------- crop frequency range --------
mask   = f >= opts.FRange(1) & f <= opts.FRange(2);
f      = f(mask);
psd    = psd(:, mask);
pwr    = pwr(:, mask);
psd_seg= psd_seg(:, mask, :);

% -------- aperiodic 1/f fit and removal --------
psd_norm = []; pwr_norm = []; psd_fit = [];
if opts.AperiodicEnable
    fmin_fit = max(1, opts.FRange(1)); if ~isempty(opts.AperiodicFminFit), fmin_fit = max(opts.AperiodicFminFit, opts.FRange(1)); end
    fmax_fit = opts.FRange(2);          if ~isempty(opts.AperiodicFmaxFit), fmax_fit = min(opts.AperiodicFmaxFit, opts.FRange(2)); end
    fitMask  = (f >= max(fmin_fit, eps)) & (f <= fmax_fit);
    if ~isempty(opts.AperiodicExcludeBands)
        ex = false(size(f));
        for ii = 1:size(opts.AperiodicExcludeBands,1)
            ex = ex | (f >= opts.AperiodicExcludeBands(ii,1) & f <= opts.AperiodicExcludeBands(ii,2));
        end
        fitMask = fitMask & ~ex;
    end
    fitMask = fitMask(:).';
    logf_fit = log10(f(fitMask));
    psd_fit = nan(size(psd), 'like', psd);
    X = [ones(numel(logf_fit),1) logf_fit(:)];

    for c = 1:nChan
        y = psd(c, fitMask); y(y<=0) = eps; logy = log10(y(:));
        if strcmpi(opts.AperiodicMethod,'robust')
            b = X \ logy;
            for iter = 1:5
                r = logy - X*b;
                s = 1.4826*mad(r,1); if s==0, s = eps; end
                w = huberWeights(r/s); W = diag(w);
                b = (X'*W*X) \ (X'*W*logy);
            end
        else
            b = X \ logy;
        end
        logSfit = b(1) + b(2)*log10(f);
        psd_fit(c,:) = 10.^logSfit;
    end

    psd_osc = psd - psd_fit; 
    psd_osc(psd_osc < 0) = 0;
    psd_norm = psd_osc;
    if numel(f) > 1, df_local = mean(diff(f)); else, df_local = df; end
    pwr_norm = psd_norm .* df_local;
end

% -------- coherence using the global window only --------
if opts.CompCoherence
    winSize = max(1, round(fs*opts.WinLength));
    step    = max(1, winSize - floor(opts.Overlap*winSize));
    segIdx  = 1:step:(nSamples - winSize + 1); if isempty(segIdx), segIdx = 1; end
    nSeg    = numel(segIdx);
    nfftC   = ternary(isempty(opts.NFFT), winSize, opts.NFFT);
    Kc      = floor(nfftC/2) + 1;

    coh_seg  = nan(nChan, nChan, Kc, nSeg);
    pcoh_seg = nan(nChan, nChan, Kc, nSeg);
    gpdc_seg = nan(nChan, nChan, Kc, nSeg);
    Su = eye(nChan);

    % one shot to get frequency axis
    idx0 = segIdx(1):(segIdx(1)+winSize-1);
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, cf0] = fdMVAR_5order(data(:,idx0), Su, nfftC, fs);
    coh_f = cf0(:).';

    % progress
    progress_init('Coherence', nSeg, usePB, usePar);

    if usePar
        q = parallel.pool.DataQueue; afterEach(q,@(~)progress_update()); %#ok<NASGU>
        parfor iSeg = 1:nSeg
            idx = segIdx(iSeg):(segIdx(iSeg)+winSize-1);
            xw  = data(:, idx);
            [~, ~, ~, gpdc, ~, coh, pcoh] = fdMVAR_5order(xw, Su, nfftC, fs);
            coh_seg(:,:,:,iSeg)  = abs(coh);
            pcoh_seg(:,:,:,iSeg) = abs(pcoh);
            gpdc_seg(:,:,:,iSeg) = abs(gpdc);
            send(q, 1);
        end
        progress_finish();
    else
        for iSeg = 1:nSeg
            idx = segIdx(iSeg):(segIdx(iSeg)+winSize-1);
            xw  = data(:, idx);
            [~, ~, ~, gpdc, ~, coh, pcoh] = fdMVAR_5order(xw, Su, nfftC, fs);
            coh_seg(:,:,:,iSeg)  = abs(coh);
            pcoh_seg(:,:,:,iSeg) = abs(pcoh);
            gpdc_seg(:,:,:,iSeg) = abs(gpdc);
            progress_step(iSeg, nSeg, usePB);
        end
        progress_finish();
    end

    if ~isempty(coh_f)
        cmask   = coh_f >= opts.FRange(1) & coh_f <= opts.FRange(2);
        coh_f   = coh_f(cmask);
        coh_seg = coh_seg(:,:,cmask,:);
        pcoh_seg= pcoh_seg(:,:,cmask,:);
        gpdc_seg= gpdc_seg(:,:,cmask,:);
    end
end

% -------- visualization --------
if opts.PlotPSD
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

if opts.PlotBands && ~isempty(psd_seg)
    figure('color','w')
    subplot(4,1,1); hold on
    idx = f <= 3;
    data_ts = squeeze(median(median(psd_seg(:,idx,:),2,'omitnan'),1,'omitnan'));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Delta'); ylabel('\muV^2/Hz'); box on

    subplot(4,1,2); hold on
    idx = f > 3 & f < 8;
    data_ts = squeeze(median(median(psd_seg(:,idx,:),2,'omitnan'),1,'omitnan'));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Theta'); ylabel('\muV^2/Hz'); box on

    subplot(4,1,3); hold on
    idx = f >= 8 & f < 13;
    data_ts = squeeze(median(median(psd_seg(:,idx,:),2,'omitnan'),1,'omitnan'));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Alpha'); ylabel('\muV^2/Hz'); box on

    subplot(4,1,4); hold on
    idx = f >= 13 & f < 30;
    data_ts = squeeze(median(median(psd_seg(:,idx,:),2,'omitnan'),1,'omitnan'));
    plot(t_seg, data_ts, '.'); plot(t_seg, smooth(data_ts, 30), 'LineWidth',2);
    title('Beta'); ylabel('\muV^2/Hz'); xlabel('Time (s)'); box on

    set(findall(gcf,'type','axes'),'FontSize',10,'FontWeight','normal');
end
end

%% Helper functions

function A = combine_segments(X, dim, how)
switch lower(how)
    case 'median', A = median(X, dim, 'omitnan');
    case 'mean',   A = mean(X, dim, 'omitnan');
    otherwise, error('Combine must be ''median'' or ''mean''.');
end
end

function [Sxx_seg, Pbin_seg, fvec, t_ends, df] = estimate_psd_fixedwin( ...
        x, fsL, winSec, method, nfftUser, overlap, nSamples, nChan, usePB, usePar, label)
% Compute PSD over fixed windows. No detrend. Parallel supported.

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
if K > 1, df = fvec(2) - fvec(1); else, df = fsL/nfft; end

% prealloc
Sxx_seg = nan(nChan, K, nSeg, 'like', x);
Pbin_seg= nan(nChan, K, nSeg, 'like', x);
t_ends  = nan(nSeg,1);

% tapers
if strcmpi(method,'multitaper')
    NW = 2.5;                          % default NW
    Kt = max(1, floor(2*NW - 1));      % default number of tapers
    [tapers, ~] = dpss(winSize, NW, Kt);   % [winSize x Kt]
else
    tapers = hamming(winSize);             % [winSize x 1]
end

% progress
progress_init(sprintf('PSD %s', lower(method)), nSeg, usePB, usePar);

if usePar
    q = parallel.pool.DataQueue; afterEach(q,@(~)progress_update()); %#ok<NASGU>
    parfor iSeg = 1:nSeg
        [Sxx_i, Pbin_i, tend_i] = psd_one_segment(x, segIdx(iSeg), winSize, fsL, nfft, K, method, tapers);
        Sxx_seg(:,:,iSeg) = Sxx_i;
        Pbin_seg(:,:,iSeg)= Pbin_i;
        t_ends(iSeg) = tend_i;
        send(q, 1);
    end
    progress_finish();
else
    for iSeg = 1:nSeg
        [Sxx_i, Pbin_i, tend_i] = psd_one_segment(x, segIdx(iSeg), winSize, fsL, nfft, K, method, tapers);
        Sxx_seg(:,:,iSeg) = Sxx_i;
        Pbin_seg(:,:,iSeg)= Pbin_i;
        t_ends(iSeg) = tend_i;
        progress_step(iSeg, nSeg, usePB);
    end
    progress_finish();
end
end

function [Sxx, Pbin, t_end] = psd_one_segment(x, startIdx, winSize, fsL, nfft, K, method, tapers)
idx   = startIdx:(startIdx+winSize-1);
t_end = idx(end)/fsL;
xseg  = x(:, idx);  % no detrend

if strcmpi(method,'multitaper')
    Sacc = zeros(size(xseg,1), K, 'like', xseg);
    for kT = 1:size(tapers,2)
        w   = tapers(:,kT).';
        W2  = sum(w.^2);
        Xw  = xseg .* w;
        Y   = fft(Xw, nfft, 2);
        Yp  = Y(:, 1:K);
        Sxx = (1/(fsL*W2)) * abs(Yp).^2;
        if rem(nfft,2) == 0, pos = 2:K-1; else, pos = 2:K; end
        Sxx(:,pos) = 2*Sxx(:,pos);
        Sacc = Sacc + Sxx;
    end
    Sxx = Sacc / size(tapers,2);
else
    w   = tapers(:).';
    W2  = sum(w.^2);
    Xw  = xseg .* w;
    Y   = fft(Xw, nfft, 2);
    Yp  = Y(:, 1:K);
    Sxx = (1/(fsL*W2)) * abs(Yp).^2;
    if rem(nfft,2) == 0, pos = 2:K-1; else, pos = 2:K; end
    Sxx(:,pos) = 2*Sxx(:,pos);
end

Pbin = Sxx * (fsL/nfft);
end

function [d, idx] = findClosest(vec, val)
[~, idx] = min(abs(vec - val)); d = vec(idx);
end

function w = huberWeights(r)
c = 1.345; a = abs(r); w = ones(size(r)); w(a > c) = c ./ a(a > c);
end

function tf = hasParallel()
tf = ~isempty(ver('parallel'));
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end

% ------------- progress utilities -------------
function progress_init(label, nTotal, usePB, usePar)
persistent p_total p_count p_label p_usePB p_usePar
p_total = max(1,nTotal); p_count = 0; p_label = label;
p_usePB = usePB; p_usePar = usePar; %#ok<NASGU>
if p_usePB, progressbar(0); end
if p_usePar, fprintf('%s: %3.0f%%', p_label, 0); end
end

function progress_update()
persistent p_total p_count p_label p_usePB p_usePar
p_count = p_count + 1;
frac = min(p_count/max(1,p_total), 1);
if p_usePB, progressbar(frac); end
if p_usePar, fprintf('\r%s: %3.0f%%', p_label, round(100*frac)); end
end

function progress_step(i, n, usePB)
if usePB, progressbar(i/n); end
end

function progress_finish()
persistent p_usePB p_usePar p_label
if ~isempty(p_usePB) && p_usePB, progressbar(1); end
if ~isempty(p_usePar) && p_usePar, fprintf('\r%s: 100%%\n', p_label); end
end
