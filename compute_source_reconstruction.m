
function source_voxel_data = compute_source_reconstruction(data, fs)

BSL = eeg_regepochs(BSL, 2, [0 2]);  % epoch in 2-s windows
leadfield   = EEG.dipfit.sourcemodel;
sourceModel = '/Users/cedriccannard/Documents/MATLAB/eeglab/functions/supportfiles/head_modelColin27_5003_Standard-10-5-Cap339.mat';
sourceModel2MNI = [0 -24 -45 0 0 -1.5708 1000 1000 1000];
nPCA = 3;
sourcemodelatlas = 'LORETA-Talairach-BAs';
headmodel = '/Users/cedriccannard/Documents/MATLAB/eeglab/plugins/dipfit/standard_BEM/standard_vol.mat';
model  = 'LCMV';  % 'LCMV' (default), 'eLoreta', 'fieldtrip' (estimation)
modelparams = {[0.0500]};


% align with MNI coordinates
tf = traditionaldipfit(g.sourcemodel2mni);
hm = load(sourceModel); % Cortex mesh or volume
pos      = tf*[hm.cortex.vertices ones(size(hm.cortex.vertices,1),1)]';
pos      = pos';
cortex.Vertices = pos(:,1:3);
cortex.Faces = hm.cortex.faces;

% make Alejandro Atlas definition compatible with Brainstrom one
nROI = length(hm.atlas.label);
cortex.Atlas(1).Name = sourcemodelatlas;
for iROI = 1:nROI
    indVertices = find(hm.atlas.colorTable == iROI);
    cortex.Atlas(1).Scouts(iROI).Label    = hm.atlas.label{iROI};
    cortex.Atlas(1).Scouts(iROI).Vertices = indVertices;
end

% Select Atlas
found = false;
for iAtlas = 1:length(cortex.Atlas)
    if strcmpi(cortex.Atlas(iAtlas).Name, sourcemodelatlas)
        cortex.Atlas = cortex.Atlas(iAtlas);
        found = true;
        break
    end
end
if ~found
    error('Atlas not found');
end

% leadfield matrix (Brainstorm or Fieldtrip)
% fieldtrip
oldLeadfield = leadfield;
leadfield.gain = reshape( [ leadfield.leadfield{:} ], [length(leadfield.label) 3 length(leadfield.leadfield)]);
leadfield.gain = permute(leadfield.gain, [1 3 2]);
leadfield = leadfield.gain;

% check number of voxels in source model match leadfield
nvox = size(cortex.Vertices, 1);
nvox2 = size(leadfield,2);
if ~isequal(nvox, nvox2)
    error('There must be the same number of vertices/voxels in the leadfield and source model');
end

% % Select channels
% chansel = 1:EEG.nbchan;
% if ~isequal(size(leadfield,1), length(chansel))
%     error('There must be the same number of channels in the leadfield and in the list of selected channels');
% end

% Frequency resolution and vector
fres = EEG.pnts/2;
frqs = sfreqs(fres, EEG.srate);

% Common average reference transform
% disp("performing common average reference to EEG data and leadfield.")
% nbchan = length(g.chansel);
% H = eye(nbchan) - ones(nbchan) ./ nbchan;
% tmpData = reshape(H*EEG.data(chansel, :), nbchan, EEG.pnts, EEG.trials); % apply to data
% leadfield = reshape(H*leadfield(:, :), nbchan, nvox, 3); % apply to leadfield

tmpData = EEG.data; % no re-reference

% Source reconstruction
if strcmpi(model, 'eLoreta')
    P_eloreta = mkfilt_eloreta_v2(leadfield, modelparams); % eLORETA inverse projection kernel
    source_voxel_data = reshape(tmpData(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3); % project to source space

elseif strcmpi(model, 'LCMV')
    C = cov(tmpData(:, :)');
    alpha = modelparams * trace(C) / length(C);
    Cr = C + alpha*eye(nbchan);
    [~, P_eloreta] = lcmv(Cr, leadfield, struct('alpha', 0, 'onedim', 0));
    source_voxel_data = reshape(tmpData(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
    source_voxel_data = 10^3*source_voxel_data; % the units are nA*m
end

%  % Transform the data back to continuous so we can get an estimate for each sample
% else
%     EEG2 = EEG;
%     EEG2.data = EEG2.data(:,:);
%     EEG2.pnts = size(EEG2.data,2);
%     EEG2.trials = 1;
%     EEG2 = eeg_checkset(EEG2);
%     dataPre = eeglab2fieldtrip(EEG2, 'preprocessing', 'dipfit');  
% 
%     % prepare data
%     cfg = [];
%     cfg.channel = {'all', '-EOG1'};
%     cfg.reref = 'yes';
%     cfg.refchannel = {'all', '-EOG1'};
%     dataPre = ft_preprocessing(cfg, dataPre);
% 
%     % load head model and prepare leadfield matrix
%     vol = load('-mat', g.headmodel);
% 
%     % source reconstruction
%     cfg             = [];
%     if lower(g.model(1)) == 'e'
%         cfg.method      = 'eLoreta';
%     else
%         cfg.method      = 'lcmv';
%     end
%     try
%         cfg.(g.sourcemethod) = struct(g.modelparams{:});
%     catch, end
%     cfg.sourcemodel = oldLeadfield;
%     cfg.headmodel   = vol.vol;
%     cfg.keeptrials  = 'yes';
%     source          = ft_sourceanalysis(cfg, dataPre);  % compute the source
% 
%     % reformat for ROI analysis below
%     source_voxel_data = reshape([ source.avg.mom{:} ], 3, size(source.avg.mom{1},2), length(source.avg.mom));
%     source_voxel_data = permute(source_voxel_data, [2 3 1]);
% end


% ROIs
nROI  = length(cortex.Atlas.Scouts);
labels = {cortex.Atlas.Scouts.Label};   % ROI labels

% keep only the strongest components (defined by nPCA) for each ROI
nfreq = fres + 1;
tmpData = reshape(source_voxel_data, EEG.pnts, EEG.trials*size(source_voxel_data,2)*size(source_voxel_data,3));
source_roi_data = [];
data_pnts = EEG.pnts;

% zero padding if necessary
if freqresolution ~= 0 % user input
    required_zeros = freqresolution - fres;
    if required_zeros < 0, error('Desired frequency resolution cannot be lower than the actual resolution of the signal.'); end
    pad = zeros(required_zeros, size(tmpData, 2));
    tmpData = cat(1,pad,tmpData,pad);
    frqs = sfreqs(freqresolution, fs);
    data_pnts = size(tmpData, 1);
    nfreq = data_pnts/2 + 1;
end

% Compute Welch power
source_roi_power = zeros(nfreq, nROI);
% size tmpData: ans = 512     4487691
% data_pnts: 512
% floor(data_pnts/2): 256
% data_pnts = 512
% EEG.srate: 256
[tmpWelch, ftmp] = pwelch(tmpData, data_pnts, floor(data_pnts/2), data_pnts, EEG.srate); % ftmp should be equal frqs 
tmpWelch = reshape(tmpWelch, size(tmpWelch,1), EEG.trials, size(source_voxel_data,2), size(source_voxel_data,3));
tmpWelch = squeeze(mean(tmpWelch,2)); % remove trials size freqs x voxels x 3
tmpWelch = squeeze(mean(tmpWelch,3)); % remove 3rd dim size freqs x voxels
source_voxel_power = tmpWelch;

% fooof settings
if strcmpi(g.fooof, 'on')
    f_range = g.fooof_frange; % freq range where 1/f should be fitted 
    settings = struct(); % use defaults
    slope = zeros(1, nROI);
    PS_corrected = zeros(size(frqs, 1), size(frqs, 2), nROI);
end
    
source_roi_data = zeros(size(source_voxel_data,1), g.nPCA*nROI);
fprintf('Computing ROI activity:');
for iROI = 1:nROI
    if mod(iROI, 5) == 0
        fprintf('.');
    end
    ind_roi = cortex.Atlas.Scouts(iROI).Vertices;
    [~, source_roi_power_norm(iROI)] = roi_getpower(source_voxel_data, ind_roi); 
    source_roi_power(:,iROI) = mean(tmpWelch(:, ind_roi),2); % shape: (101, nROI)

    [source_roi_data_tmp, nPCA(iROI)] = roi_getact(source_voxel_data, ind_roi, g.nPCA);
    source_roi_data(:, (iROI-1)*g.nPCA+1:iROI*g.nPCA) = source_roi_data_tmp;
    if strcmpi(g.fooof, 'on')
        ps1 = source_roi_power(:,iROI);
        fooof_result = fooof(frqs, ps1, f_range, settings, true);
        
        offset = fooof_result.aperiodic_params(1);
        slope(iROI) = fooof_result.aperiodic_params(2);
        y = (-slope(iROI) .* log10(frqs)) + offset;
        PS_corrected(:,:,iROI) = 10*log10(ps1)-10*y;
    end
end
fprintf('\n');

% version with nPCA components
source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);

