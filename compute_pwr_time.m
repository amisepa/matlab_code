function [psd, pwr_t, f] = compute_pwr_time(data, fs, freqHz, varargin)
% COMPUTE_PWR_TIME  PSD over time using Welch, preserving frequency bins.
%
% Inputs:
%   data   : [nChan x nSamples] continuous signal
%   fs     : sampling rate (Hz)
%   freqHz : [fmin fmax] overall frequency range to return (Hz)
%
% Outputs:
%   psd   : [nFreq x nSeg] PSD over time within freqHz (uV^2/Hz or dB)
%   pwr_t : [nSeg x 1] time vector (s), segment end times
%   f     : [nFreq x 1] frequency vector (Hz), limited to freqHz
%
% Name-value options:
%   'Chan'        : [] = all channels (20% trimmed mean),
%                   numeric indices, or logical mask length nChan
%   'DB'          : true for 10*log10 (default false)
%   'Win'         : window length in seconds (default 2)
%   'Plot'        : true/false plot (default true)
%   'Smooth'      : smoothing window (points) for plotting only (default 15)
%   'UseParallel' : true/false use parfor across segments (default true)
%   'Progress'    : true/false progress display (default true)
%
% Cedric Cannard, Jan 2026

p = inputParser;
p.addParameter('Chan', [], @(x) isempty(x) || isnumeric(x) || islogical(x));
p.addParameter('DB', false, @(x) islogical(x) || ismember(x,[0 1]));
p.addParameter('Win', 2, @(x) isscalar(x) && x > 0);
p.addParameter('Plot', true, @(x) islogical(x) || ismember(x,[0 1]));
p.addParameter('Smooth', 15, @(x) isscalar(x) && x >= 1);
p.addParameter('UseParallel', true, @(x) islogical(x) || ismember(x,[0 1]));
p.addParameter('Progress', true, @(x) islogical(x) || ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

if nargin < 3 || isempty(freqHz) || ~isnumeric(freqHz) || numel(freqHz) ~= 2 || any(~isfinite(freqHz))
    error('freqHz must be a 2-element numeric vector: [fmin fmax].');
end
freqHz = sort(freqHz(:).');
if freqHz(1) < 0 || freqHz(2) <= freqHz(1)
    error('freqHz must satisfy 0 <= fmin < fmax.');
end

if size(data,1) > size(data,2)
    data = data.';
end
[nChan, nSamples] = size(data);

if ~isempty(opt.Chan)

    if islogical(opt.Chan)
        if ~isvector(opt.Chan) || numel(opt.Chan) ~= nChan
            error('Chan logical mask must be a vector of length %d (nChan).', nChan);
        end
        opt.Chan = find(opt.Chan(:)).';
    else
        if ~isvector(opt.Chan) || any(~isfinite(opt.Chan)) || any(opt.Chan ~= round(opt.Chan)) || any(opt.Chan < 1)
            error('Chan must be [] or a vector of positive integer indices, or a logical mask.');
        end
        if any(opt.Chan > nChan)
            error('Chan indices exceed number of channels (%d).', nChan);
        end
        opt.Chan = opt.Chan(:).';
    end

    if isempty(opt.Chan)
        error('Chan mask/indices selected 0 channels.');
    end
end

if opt.Progress
    if isempty(opt.Chan)
        fprintf('Channels: all (%d), using 20%% trimmed mean\n', nChan);
    elseif isscalar(opt.Chan)
        fprintf('Channel: %d\n', opt.Chan);
    else
        fprintf('Channels: %d selected, averaging\n', numel(opt.Chan));
    end
end

winSamp  = min(round(opt.Win * fs), nSamples);
noverlap = floor(0.5 * winSamp);
nfft     = winSamp;

step   = winSamp - noverlap;
segBeg = 1:step:(nSamples - winSamp + 1);
nSeg   = numel(segBeg);

pwr_t = nan(nSeg,1);
w = hamming(winSamp);

[~, f_all] = pwelch(data(1,1:winSamp).', w, noverlap, nfft, fs);
maskF = f_all >= freqHz(1) & f_all < freqHz(2);
f = f_all(maskF);
nFreq = numel(f);
if nFreq < 2
    error('freqHz selects too few bins. Increase range or window length.');
end

psd_ch = nan(nChan, nFreq, nSeg);

usePar = opt.UseParallel && ~isempty(ver('parallel'));
if usePar && isempty(gcp('nocreate'))
    parpool;
end
usePar = usePar && ~isempty(gcp('nocreate'));

if opt.Progress
    fprintf('Computing PSD (%d segments)...\n', nSeg);
end

if usePar
    dq = parallel.pool.DataQueue;
    nDone = 0;
    if opt.Progress
        afterEach(dq, @updateProgress);
    end

    parfor iSeg = 1:nSeg
        idx = segBeg(iSeg):(segBeg(iSeg)+winSamp-1);
        pwr_t(iSeg) = idx(end) / fs;

        Pxx = pwelch(data(:,idx).', w, noverlap, nfft, fs); % [nFreqAll x nChan]
        psd_ch(:,:,iSeg) = Pxx(maskF,:).';                 % [nChan x nFreq x nSeg]

        if opt.Progress
            send(dq, iSeg);
        end
    end
else
    for iSeg = 1:nSeg
        idx = segBeg(iSeg):(segBeg(iSeg)+winSamp-1);
        pwr_t(iSeg) = idx(end) / fs;

        Pxx = pwelch(data(:,idx).', w, noverlap, nfft, fs);
        psd_ch(:,:,iSeg) = Pxx(maskF,:).';

        if opt.Progress && mod(iSeg, max(1, round(nSeg/10))) == 0
            fprintf('  %3.0f%% complete\n', 100*iSeg/nSeg);
        end
    end
end

if isempty(opt.Chan)
    psd = squeeze(trimmean(psd_ch, 20, 1, 'round')); % [nFreq x nSeg]
    chan_label = 'all chans (20% trimmean)';
else
    psd = squeeze(mean(psd_ch(opt.Chan,:,:), 1, 'omitnan')); % [nFreq x nSeg]
    if isscalar(opt.Chan)
        chan_label = sprintf('chan %d', opt.Chan);
    else
        chan_label = sprintf('%d chans (mean)', numel(opt.Chan));
    end
end

if opt.DB
    psd = 10 * log10(max(psd, eps));
    zlab = 'PSD (dB)';
else
    zlab = 'PSD (\muV^2/Hz)';
end

if opt.Plot
    df = mean(diff(f));
    pwr = sum(max(psd, 0), 1).' * df; % only for plotting summary

    t_min = pwr_t / 60;
    figure('Color','w');
    plot(t_min, pwr, '.'); hold on

    wsm = min(opt.Smooth, numel(pwr));
    if wsm >= 2
        pwr_sm = movmean(pwr, wsm, 'omitnan');
        k = floor(wsm/2);
        ii = (1+k):(numel(pwr)-k);
        plot(t_min(ii), pwr_sm(ii), 'LineWidth', 1.5);
        legend({'raw','smoothed'}, 'Location','best');
    end
    xlabel('Time (min)');
    ylabel(sprintf('Total power in %.1f-%.1f Hz (a.u.)', freqHz(1), freqHz(2)));
    title(sprintf('Welch PSD over time | %s | %s', chan_label, zlab));
    grid on; box on
end

if opt.Progress
    fprintf('Done.\n');
end

    function updateProgress(~)
        nDone = nDone + 1;
        if mod(nDone, max(1, round(nSeg/10))) == 0 || nDone == nSeg
            fprintf('  %3.0f%% complete\n', 100*nDone/nSeg);
        end
    end
end
