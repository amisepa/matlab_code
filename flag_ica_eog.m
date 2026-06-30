function [isOcular, scores, details] = flag_ica_eog(EEG, EOG, varargin)
% FLAG_ICA_EOG  Flag ocular ICA components using correlation or regression.
%
% Purpose
%   Identify ICA components that track ocular activity when ICLabel is unreliable
%   or you prefer a signal-based detector. Works well with small EEG montages.
%
% Inputs
%   EEG  : EEGLAB EEG struct with ICA done.
%          Uses EEG.icaact if present, otherwise computes from icaweights,
%          icasphere, and EEG.data(EEG.icachansind,:).
%   EOG  : Either an EEGLAB dataset with 1 or 2 EOG channels, or a numeric matrix
%          of size [nSamples x 1 or 2]. If one EOG channel is provided it is
%          duplicated internally so both sum and difference references exist.
%
% Name-Value pairs
%   'Method'      : 'correlation' | 'regression' | 'both'  (default 'correlation')
%   'Band'        : [] or [lo hi] in Hz for optional bandpass of blinks and saccades.
%                   Default [] (no filtering). Typical choices:
%                   [0.03 6] very slow large blinks,
%                   [0.05 8] reliable default for blinks,
%                   [0.1 12] blinks plus saccades,
%                   [1 15] faster eye activity when drift is strong.
%   'LagSec'      : Maximum absolute lag for cross correlation in seconds.
%                   Default 0.8.
%                   0.3  tighter match when timing is well aligned;
%                   1.2  more tolerant when EOG and EEG are slightly misaligned.
%   'ZscoreMAD'   : true uses median and MAD for robust z-scoring. false uses mean and SD.
%                   Default true.
%   'PMethod'     : 'fisher' | 'perm' | 'none'. P-value method for correlation metrics.
%                   'fisher' uses Fisher z with Bonferroni over lags and EOG refs. Fast.
%                   'perm' uses circular shift permutations that preserve EOG autocorrelation.
%                   'none' skips p-values and uses only effect-size floors.
%                   Default 'fisher'.
%   'MinAbsR'     : Effect-size floor for correlation pathway. Default 0.05.
%                   0.10 stricter, fewer components flagged;
%                   0.02 more sensitive, may catch subtle saccades.
%   'MinR2'       : Effect-size floor for regression pathway. Default 0.05.
%                   0.10 stricter, flags only strong ocular ICs;
%                   0.02 more sensitive, may include weaker matches.
%   'Nperm'       : Number of permutations if PMethod is 'perm'. Default 800.
%                   100 faster but less stable p values;
%                   800 slower but more stable p values.
%   'ReportTop'   : If > 0, prints and plots the top K ICs by score. Default 1.
%   'OpenProp'    : If true and ReportTop > 0, calls pop_prop on the top K ICs. Default false.
%   'UseParallel' : true to parallelize across ICs with parfor. Default false.
%   'ThreshP'     : FDR alpha across components. Default 0.05.
%
% Outputs
%   isOcular : Logical [nComp x 1]. True if component is flagged for removal.
%   scores   : Table with one row per component. Columns:
%              IC, score, rho_spearman, rho_xcorr, r2_reg, p_perm_r2,
%              p_perm_corr_min, p_final, isOcular
%   details  : Struct with arrays and settings for auditing.

% ---------------- options
p = inputParser;
p.addParameter('Method','correlation');
p.addParameter('Band',[]);
p.addParameter('LagSec',0.8);
p.addParameter('ZscoreMAD',true);
p.addParameter('MinAbsR',0.05);
p.addParameter('MinR2',0.05);
p.addParameter('PMethod','fisher');
p.addParameter('Nperm',800);
p.addParameter('ReportTop',1);
p.addParameter('OpenProp',false);
p.addParameter('UseParallel',false);
p.addParameter('ThreshP',0.05);
p.parse(varargin{:});
opt   = p.Results;
method = lower(opt.Method);
pmeth  = lower(opt.PMethod);

% pretty log suffix
if opt.UseParallel, logSuffix = ', parallel'; else, logSuffix = ''; end

% ---------------- IC activations, time in rows
assert(isfield(EEG,'icaweights') && ~isempty(EEG.icaweights),'EEG has no ICA weights.');
if ~isfield(EEG,'icaact') || isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind,:);
end
X = double(EEG.icaact)';     % [nSamp x nComp]
[nSamp, nComp] = size(X);
fs = EEG.srate;

% ---------------- EOG to [nSamp x nEOG]
if isstruct(EOG) && isfield(EOG,'data')
    E = double(EOG.data)';   % [nSamp x nEOG]
    fsE = EOG.srate;
elseif isnumeric(EOG)
    E = double(EOG);         % assume [nSamp x nEOG]
    fsE = fs;
else
    error('EOG must be EEGLAB struct or numeric matrix.');
end
if isempty(E) || size(E,2) == 0, error('EOG has zero channels.'); end
if size(E,2) == 1, E = [E E]; end
if fsE ~= fs, E = resample(E, fs, fsE); end

% align by truncation
L = min(size(E,1), nSamp);
X = X(1:L,:); E = E(1:L,:); nSamp = L;

% EOG refs: E1, E2, sum, diff
E1 = E(:,1); E2 = E(:,2);
EOGX = [E1, E2, (E1 + E2), (E1 - E2)];   % [nSamp x 4]
nRef = size(EOGX,2);

% ---------------- filtering and scaling
EOGXf = EOGX; Xf = X;
ordMax  = 4;
ordSafe = min(ordMax, floor((nSamp - 1)/6));
doFilt  = ~isempty(opt.Band) && ordSafe >= 1;
if doFilt
    [b, a] = butter(ordSafe, opt.Band/(fs/2), 'bandpass');
    for k = 1:nRef,  EOGXf(:,k) = filtfilt(b, a, detrend(EOGX(:,k),'constant')); end
    for c = 1:nComp, Xf(:,c)    = filtfilt(b, a, detrend(X(:,c),   'constant')); end
else
    for k = 1:nRef,  EOGXf(:,k) = detrend(EOGX(:,k),'constant'); end
    for c = 1:nComp, Xf(:,c)    = detrend(X(:,c),   'constant'); end
end
if opt.ZscoreMAD
    EOGXf = zscore_mad(EOGXf);
    Xf    = zscore_mad(Xf);
else
    EOGXf = zscore(EOGXf);
    Xf    = zscore(Xf);
end

% xcorr lag and indices
lagMax = max(0, min( round(opt.LagSec*fs), floor(nSamp/4) ));
lags   = -lagMax:lagMax;
nLag   = numel(lags);
Lcorr  = 2*nSamp - 1;
center = nSamp;
win    = (-lagMax:lagMax);
baseIdx = center + win;

% ---------------- allocate
rho_spear   = nan(nComp, nRef);
rho_xcorr   = nan(nComp, nRef);
p_perm_corr = nan(nComp, nRef);
r2_reg      = nan(nComp, 1);
p_perm_r2   = nan(nComp, 1);
r2_null_p95 = nan(nComp, 1);

% optional pool
if opt.UseParallel
    try
        pool = gcp('nocreate');
        if isempty(pool), parpool; end
    catch
        warning('UseParallel requested but parpool could not start. Falling back to serial.');
        opt.UseParallel = false; logSuffix = '';
    end
end

% permutation shifts reused across ICs
if strcmp(pmeth,'perm') && opt.Nperm > 0
    perm_shifts = randi(nSamp, opt.Nperm, 1) - 1;  % 0..nSamp-1
else
    perm_shifts = [];
end

% ---------------- correlation pathway  (fast perms via shifted xcorr window)
doCorr = any(strcmp(method, {'correlation','both'}));
if doCorr
    fprintf('Correlation path: %d ICs, %d EOG refs, %d lags%s\n', ...
        nComp, nRef, nLag, logSuffix);

    if ~isempty(perm_shifts)
        idxMatPerm = mod(baseIdx.' + perm_shifts.' - 1, Lcorr) + 1; % [nLag x Nperm]
    else
        idxMatPerm = [];
    end

    P = startProgress('Correlation', nComp, opt.UseParallel);

    if opt.UseParallel
        parfor c = 1:nComp
            xc = Xf(:,c);
            loc_s = zeros(1,nRef);
            loc_x = zeros(1,nRef);
            loc_p = nan(1,nRef);
            for j = 1:nRef
                ej = EOGXf(:,j);
                loc_s(j) = abs(corr(xc, ej, 'Type','Spearman', 'Rows','complete'));
                cf = xcorr(xc, ej, 'coeff');
                rObs = max(abs(cf(baseIdx)));
                loc_x(j) = rObs;
                if ~isempty(idxMatPerm)
                    nulls = max(abs(cf(idxMatPerm)), [], 1).';
                    loc_p(j) = max(1, sum(nulls >= rObs)) / numel(perm_shifts);
                end
            end
            rho_spear(c,:) = loc_s;
            rho_xcorr(c,:) = loc_x;
            p_perm_corr(c,:) = loc_p;
            if ~isempty(P.dq), send(P.dq,1); else, P.update(); end
        end
    else
        for c = 1:nComp
            xc = Xf(:,c);
            for j = 1:nRef
                ej = EOGXf(:,j);
                rho_spear(c,j) = abs(corr(xc, ej, 'Type','Spearman', 'Rows','complete'));
                cf = xcorr(xc, ej, 'coeff');
                rObs = max(abs(cf(baseIdx)));
                rho_xcorr(c,j) = rObs;
                if ~isempty(idxMatPerm)
                    nulls = max(abs(cf(idxMatPerm)), [], 1).';
                    p_perm_corr(c,j) = max(1, sum(nulls >= rObs)) / numel(perm_shifts);
                end
            end
            P.update();
        end
    end
    P.close();
end

% ---------------- regression pathway  (orthogonalized refs + fast perms)
doReg = any(strcmp(method, {'regression','both'}));
Ereg = []; r = 0; Ginv = []; A = [];
if doReg
    fprintf('Regression path: %d ICs, %d EOG refs%s\n', nComp, nRef, logSuffix);

    % orthogonalize EOG refs and keep first 2 axes
    [Q,~] = qr(zscore(EOGXf,0,1), 0);
    r = min(2, size(Q,2));
    Ereg = Q(:,1:r);

    % design and Gram
    A    = [Ereg, ones(nSamp,1)];
    G    = A.' * A;
    Ginv = inv(G);

    % indices for permutations
    if ~isempty(perm_shifts)
        idxPerm = center + perm_shifts;
    else
        idxPerm = [];
    end

    P = startProgress('Regression', nComp, opt.UseParallel);

    if opt.UseParallel
        parfor c = 1:nComp
            y   = Xf(:,c);
            y0  = y - mean(y);
            SST = sum(y0.^2);
            sumy = sum(y);

            v0   = A.' * y;
            yHy  = v0.' * (Ginv * v0);
            r2o  = max(0, min(1, yHy / SST));
            pLoc = NaN; r2p95 = NaN;

            if ~isempty(idxPerm)
                cfMat = zeros(2*nSamp-1, r);
                for j = 1:r
                    cfMat(:,j) = xcorr(y, Ereg(:,j));
                end
                V = zeros(r+1, numel(idxPerm));
                for j = 1:r
                    V(j,:) = cfMat(idxPerm, j).';
                end
                V(end,:) = sumy;
                W        = Ginv * V;
                yHy_perm = sum(V .* W, 1);
                r2_perm  = max(0, min(1, yHy_perm ./ SST));
                pLoc     = max(1, sum(r2_perm >= r2o)) / numel(idxPerm);
                r2p95    = prctile(r2_perm, 95);
            end

            r2_reg(c)      = r2o;
            p_perm_r2(c)   = pLoc;
            r2_null_p95(c) = r2p95;

            if ~isempty(P.dq), send(P.dq,1); else, P.update(); end
        end
    else
        for c = 1:nComp
            y   = Xf(:,c);
            y0  = y - mean(y);
            SST = sum(y0.^2);
            sumy = sum(y);

            v0   = A.' * y;
            yHy  = v0.' * (Ginv * v0);
            r2o  = max(0, min(1, yHy / SST));
            r2_reg(c) = r2o;

            if ~isempty(idxPerm)
                cfMat = zeros(2*nSamp-1, r);
                for j = 1:r
                    cfMat(:,j) = xcorr(y, Ereg(:,j));
                end
                V = zeros(r+1, numel(idxPerm));
                for j = 1:r
                    V(j,:) = cfMat(idxPerm, j).';
                end
                V(end,:) = sumy;
                W        = Ginv * V;
                yHy_perm = sum(V .* W, 1);
                r2_perm  = max(0, min(1, yHy_perm ./ SST));
                p_perm_r2(c)   = max(1, sum(r2_perm >= r2o)) / numel(idxPerm);
                r2_null_p95(c) = prctile(r2_perm, 95);
            end

            P.update();
        end
    end
    P.close();
end

% adaptive regression floor from empirical null
MinR2_eff = opt.MinR2;
if any(isfinite(r2_null_p95))
    MinR2_eff = max(opt.MinR2, median(r2_null_p95, 'omitnan'));
end

% analytic regression p values (F-test) used when PMethod is not 'perm'
p_reg_analytic = [];
if doReg && ~strcmpi(pmeth,'perm')
    k = size(Ereg,2);              % number of regression EOG axes used
    n = nSamp;
    r2_clip = min(max(r2_reg,0), 0.999999);
    Fstat = (r2_clip./max(1,k)) ./ ((1 - r2_clip) ./ max(1, n - k - 1));
    p_reg_analytic = 1 - fcdf(Fstat, max(1,k), max(1, n - k - 1));
end

% ---------------- combine metrics
switch method
    case 'correlation'
        scoreCorr = max([rho_spear rho_xcorr], [], 2, 'omitnan');
        if strcmp(pmeth,'fisher')
            mtests = max(1, nLag * nRef);
            r = max(rho_xcorr,[],2,'omitnan');
            nz = max(4, nSamp - lagMax);
            z = atanh(max(min(r,0.999999),-0.999999)) .* sqrt(nz - 3);
            pCorr = 2 * (1 - normcdf(abs(z)));
            pCorr = min(1, pCorr * mtests);
        elseif strcmp(pmeth,'perm')
            pCorr = min(p_perm_corr, [], 2, 'omitnan');
        else
            pCorr = ones(nComp,1);
        end
        score = scoreCorr;
        p_use = pCorr;
        p_use(score < opt.MinAbsR | ~isfinite(p_use)) = 1;

    case 'regression'
        score = r2_reg;
        if strcmpi(pmeth,'perm')
            p_use = p_perm_r2;
        elseif ~isempty(p_reg_analytic)
            p_use = p_reg_analytic;
        else
            p_use = ones(nComp,1);   % effect-only gating if desired
        end
        p_use(score < MinR2_eff | ~isfinite(p_use)) = 1;

    case 'both'
        scoreCorr = max([rho_spear rho_xcorr], [], 2, 'omitnan');
        scoreReg  = r2_reg;

        % correlation p's
        if strcmp(pmeth,'fisher')
            mtests = max(1, nLag * nRef);
            r = max(rho_xcorr,[],2,'omitnan');
            nz = max(4, nSamp - lagMax);
            z = atanh(max(min(r,0.999999),-0.999999)) .* sqrt(nz - 3);
            pCorr = 2 * (1 - normcdf(abs(z)));
            pCorr = min(1, pCorr * mtests);
        elseif strcmp(pmeth,'perm')
            pCorr = min(p_perm_corr, [], 2, 'omitnan');
        else
            pCorr = ones(nComp,1);
        end

        % regression p's
        if strcmpi(pmeth,'perm')
            pReg = p_perm_r2;
        elseif ~isempty(p_reg_analytic)
            pReg = p_reg_analytic;
        else
            pReg = ones(nComp,1);
        end

        % floors
        pCorr(scoreCorr < opt.MinAbsR | ~isfinite(pCorr)) = 1;
        pReg( scoreReg  < MinR2_eff  | ~isfinite(pReg) )  = 1;

        score = max([scoreCorr scoreReg], [], 2, 'omitnan');
        p_use = min([pCorr pReg], [], 2, 'omitnan');
end

isOcular = fdr_bh_local(p_use, opt.ThreshP);

% ---------------- outputs
scores = table((1:nComp)', score, ...
    max(rho_spear,[],2,'omitnan'), max(rho_xcorr,[],2,'omitnan'), ...
    r2_reg, p_perm_r2, ...
    min(p_perm_corr,[],2,'omitnan'), p_use, isOcular, ...
    'VariableNames', {'IC','score','rho_spearman','rho_xcorr', ...
                      'r2_reg','p_perm_r2','p_perm_corr_min','p_final','isOcular'});

details = struct();
details.Method        = method;
details.PMethod       = pmeth;
details.UseParallel   = opt.UseParallel;
details.rho_spearman  = rho_spear;
details.rho_xcorr     = rho_xcorr;
details.r2_reg        = r2_reg;
details.p_perm_corr   = p_perm_corr;
details.p_perm_r2     = p_perm_r2;
details.r2_null_p95   = r2_null_p95;
details.MinR2_eff     = MinR2_eff;
details.Band          = opt.Band;
details.LagSec        = opt.LagSec;
details.MinAbsR       = opt.MinAbsR;
details.MinR2         = opt.MinR2;
details.nLag          = nLag;
details.ThreshP       = opt.ThreshP;
details.RegAxesUsed   = size(Ereg,2);

% ---------------- optional report: always show top K by score, highlight flagged
if opt.ReportTop > 0
    [~, order] = sort(scores.score, 'descend');
    K = min(opt.ReportTop, nComp);
    topIdx = order(1:K);

    fprintf('\nTop %d components by score (flagged highlighted):\n', K);
    disp(scores(topIdx, :))

    try
        figure('Name','Top ICs by score');
        stem(scores.IC(topIdx), scores.score(topIdx), 'filled'); hold on
        topFlag = isOcular(topIdx);
        if any(topFlag)
            plot(scores.IC(topIdx(topFlag)), scores.score(topIdx(topFlag)), 'o', ...
                 'MarkerSize', 8, 'LineWidth', 1.5);
            legend('Top by score','Flagged','Location','best');
        else
            legend('Top by score','Location','best');
        end
        xlabel('IC'); ylabel('score'); title(sprintf('Top %d ICs by score', K));
        grid on
    catch
        warning('Quick plot failed.');
    end

    if opt.OpenProp
        for k = 1:K
            try, pop_prop(EEG, 0, topIdx(k)); catch, warning('pop_prop failed.'); end
        end
    end
end
end

% ---------------- helpers
function Xz = zscore_mad(X)
    med  = median(X,1,'omitnan');
    madv = mad(X,1,1);
    madv(madv==0) = 1;
    Xz = (X - med) ./ madv;
end

function h = fdr_bh_local(p, q)
    ps = p(:);
    ps(~isfinite(ps)) = 1;
    [pss, idx] = sort(ps);
    m = numel(pss);
    thr = (1:m)' .* q ./ m;
    k = find(pss <= thr, 1, 'last');
    h = false(m,1);
    if ~isempty(k), h(idx(1:k)) = true; end
end

% progress utility
function P = startProgress(label, total, useParallel)
    P.total = total;
    P.count = 0;
    P.notifyEvery = max(1, round(total/20));
    P.hwb = [];
    try, P.hwb = waitbar(0, label); end
    if useParallel
        P.dq = parallel.pool.DataQueue;
        afterEach(P.dq, @(~) update());
    else
        P.dq = [];
    end
    function update()
        P.count = P.count + 1;
        if ~isempty(P.hwb) && (mod(P.count, P.notifyEvery) == 0 || P.count == P.total)
            waitbar(P.count/P.total, P.hwb, sprintf('%s %d/%d', label, P.count, P.total));
        end
        if mod(P.count, P.notifyEvery) == 0 || P.count == P.total
            fprintf('%s %d/%d (%.0f%%)\n', label, P.count, P.total, 100*P.count/P.total);
        end
    end
    P.update = @update;
    P.close  = @() ( ~isempty(P.hwb) && close(P.hwb) );
end
