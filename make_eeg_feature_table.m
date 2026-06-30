function T = make_eeg_feature_table(EEG, subject_id, chan_labels, badChan)
% MAKE_EEG_FEATURE_TABLE  Build per-channel feature table for ML
%
% Inputs
%   EEG          struct with EEG.brainbeats.features.EEG.* as you showed
%   subject_id   char or string, becomes a column
%   chan_labels  cellstr of channel names, length = nChan
%   badChan      logical vector [1 x nChan] or [nChan x 1], 1 means bad
%
% Output
%   T            table with one row per channel and numeric predictors
%
% Notes
%   1) All predictors for bad channels are set to NaN and a flag bad_chan is set.
%   2) Only scalar per-channel features are included. Full PSD is not added.
%   3) Asymmetry values are skipped by default since they are pairwise, not per-channel.

F  = EEG.brainbeats.features.EEG;
nC = numel(F.time.rms);

if nargin < 4 || isempty(badChan), badChan = false(nC,1); end
badChan = badChan(:) ~= 0;                       % column logical

% Convenience handles
Ft  = F.time;
Ff  = F.frequency;
Fq  = F.qeeg;
Fn  = F.nonlinear;

% Preallocate a struct array, one element per channel
rows(nC,1) = struct( ...
    'subject',      string(missing), ...
    'channel',      string(missing), ...
    'bad_chan',     false, ...
    ...
    'time_rms',     NaN, 'time_mode', NaN, 'time_var', NaN, ...
    'time_skew',    NaN, 'time_kurt', NaN, 'time_iqr', NaN, ...
    ...
    'band_delta',   NaN, 'band_theta', NaN, 'band_alpha', NaN, ...
    'band_beta',    NaN, 'band_gamma', NaN, ...
    ...
    'qeeg_delta_abs', NaN, 'qeeg_theta_abs', NaN, 'qeeg_alpha_abs', NaN, ...
    'qeeg_beta_abs',  NaN, 'qeeg_gamma_abs', NaN, 'qeeg_total_abs', NaN, ...
    'qeeg_delta_rel', NaN, 'qeeg_theta_rel', NaN, 'qeeg_alpha_rel', NaN, ...
    'qeeg_beta_rel',  NaN, 'qeeg_gamma_rel', NaN, ...
    'qeeg_alpha_theta', NaN, 'qeeg_theta_beta', NaN, ...
    'qeeg_alpha_beta',  NaN, 'qeeg_alpha_tplusb', NaN, ...
    'qeeg_iaf',     NaN, 'qeeg_median', NaN, 'qeeg_sef90', NaN, ...
    ...
    'nl_fd',        NaN, 'nl_fe', NaN );

% Fill rows
sid = string(subject_id);
for c = 1:nC
    rows(c).subject  = sid;
    if nargin >= 3 && ~isempty(chan_labels)
        rows(c).channel = string(chan_labels{c});
    else
        rows(c).channel = "C"+c;
    end
    rows(c).bad_chan = badChan(c);

    % If bad channel, leave NaNs and continue
    % If bad channel, overwrite every numeric predictor with NaN so table stays numeric
    if badChan(c)
        rows(c).time_rms   = NaN; rows(c).time_mode  = NaN; rows(c).time_var = NaN;
        rows(c).time_skew  = NaN; rows(c).time_kurt  = NaN; rows(c).time_iqr = NaN;
    
        rows(c).band_delta = NaN; rows(c).band_theta = NaN; rows(c).band_alpha = NaN;
        rows(c).band_beta  = NaN; rows(c).band_gamma = NaN;
    
        rows(c).qeeg_delta_abs = NaN; rows(c).qeeg_theta_abs = NaN; rows(c).qeeg_alpha_abs = NaN;
        rows(c).qeeg_beta_abs  = NaN; rows(c).qeeg_gamma_abs = NaN; rows(c).qeeg_total_abs = NaN;
    
        rows(c).qeeg_delta_rel = NaN; rows(c).qeeg_theta_rel = NaN; rows(c).qeeg_alpha_rel = NaN;
        rows(c).qeeg_beta_rel  = NaN; rows(c).qeeg_gamma_rel = NaN;
    
        rows(c).qeeg_alpha_theta = NaN; rows(c).qeeg_theta_beta = NaN;
        rows(c).qeeg_alpha_beta  = NaN; rows(c).qeeg_alpha_tplusb = NaN;
    
        rows(c).qeeg_iaf = NaN; rows(c).qeeg_median = NaN; rows(c).qeeg_sef90 = NaN;
    
        rows(c).nl_fd = NaN; rows(c).nl_fe = NaN;
        continue
    end

    % time
    rows(c).time_rms   = safeget(Ft.rms, c);
    rows(c).time_mode  = safeget(Ft.mode, c);
    rows(c).time_var   = safeget(Ft.var, c);
    rows(c).time_skew  = safeget(Ft.skewness, c);
    rows(c).time_kurt  = safeget(Ft.kurtosis, c);
    rows(c).time_iqr   = safeget(Ft.iqr, c);

    % frequency band powers
    rows(c).band_delta = safeget(Ff.delta, c);
    rows(c).band_theta = safeget(Ff.theta, c);
    rows(c).band_alpha = safeget(Ff.alpha, c);
    rows(c).band_beta  = safeget(Ff.beta,  c);
    rows(c).band_gamma = safeget(Ff.gamma, c);

    % qEEG absolute and relative
    rows(c).qeeg_delta_abs = safeget(Fq.delta_abs, c);
    rows(c).qeeg_theta_abs = safeget(Fq.theta_abs, c);
    rows(c).qeeg_alpha_abs = safeget(Fq.alpha_abs, c);
    rows(c).qeeg_beta_abs  = safeget(Fq.beta_abs,  c);
    rows(c).qeeg_gamma_abs = safeget(Fq.gamma_abs, c);
    rows(c).qeeg_total_abs = safeget(Fq.total_abs,  c);

    rows(c).qeeg_delta_rel = safeget(Fq.delta_rel, c);
    rows(c).qeeg_theta_rel = safeget(Fq.theta_rel, c);
    rows(c).qeeg_alpha_rel = safeget(Fq.alpha_rel, c);
    rows(c).qeeg_beta_rel  = safeget(Fq.beta_rel,  c);
    rows(c).qeeg_gamma_rel = safeget(Fq.gamma_rel, c);

    % qEEG ratios and indices
    rows(c).qeeg_alpha_theta = safeget(Fq.alpha_theta, c);
    rows(c).qeeg_theta_beta  = safeget(Fq.theta_beta,  c);
    rows(c).qeeg_alpha_beta  = safeget(Fq.alpha_beta,  c);
    rows(c).qeeg_alpha_tplusb= safeget(Fq.alpha_tplusb, c);
    rows(c).qeeg_iaf         = safeget(Fq.iaf, c);
    rows(c).qeeg_median      = safeget(Fq.median, c);
    rows(c).qeeg_sef90       = safeget(Fq.sef90, c);

    % nonlinear
    rows(c).nl_fd = safeget(Fn.FD, c);
    rows(c).nl_fe = safeget(Fn.FE, c);
end

% Convert to table
T = struct2table(rows);

% Make variable names valid and unique (in case labels collide)
T.Properties.VariableNames = matlab.lang.makeUniqueStrings( ...
    matlab.lang.makeValidName(T.Properties.VariableNames));

end

function x = safeget(vecOrMat, idx)
% Return element idx from a vector; if empty or too short, return NaN
try
    if isempty(vecOrMat), x = NaN; return; end
    if ismatrix(vecOrMat) && size(vecOrMat,2) == 1
        x = vecOrMat(idx,1);
    elseif ismatrix(vecOrMat) && size(vecOrMat,1) == 1
        x = vecOrMat(1,idx);
    else
        x = vecOrMat(idx); % default linear indexing for column vectors
    end
catch
    x = NaN;
end
end
