%% Detect flat lines in EEG data, excludes very small segments (< 5 samples),
% plots them in red (superimposed on raw EEG data), and print portion of data
% that were removed (seconds and ratio).
%
% Copyright (C) - Cedric Cannard, 2023

function EEG = rm_flatSegments(EEG,flatThresh,vis)

% Thresholds
if ~exist('flatThresh','var') || isempty(flatThresh)
    % flatThresh = 20*eps;
    flatThresh = 0.05;
end

% Visualization
if ~exist('vis','var') || isempty(vis)
    vis = true;
end

% copy of original dataset for visualization
if vis
    oriEEG = EEG;
end

nChan = EEG.nbchan;
mask = nan(nChan,EEG.pnts+1);
for iChan = 1:nChan
    mask(iChan,:) = [false abs(diff(EEG.data(iChan,:)))<flatThresh false]; % false 1st and last sample to preserve original size
    % mask(iChan,:) = [false isoutlier(abs(diff(EEG.data(iChan,:))),'mean') false];   % this can detect EEG artifacts
end
if nChan == 4
    mask_front = or(mask(2,:),mask(3,:));
    mask_post = or(mask(1,:),mask(4,:));
    mask = or(mask_front, mask_post);
elseif nChan == 3
    flat_mask_tmp = or(mask(1,:),mask(2,:));
    mask = or(mask(3,:),flat_mask_tmp);
elseif nChan == 2
    mask = or(mask(1,:),mask(2,:));
end
mask(end) = [];

% Remove them if present
if sum(mask) > 0
    flatSeg = reshape(find(diff([false mask false])),2,[])';
    flatSeg(:,2) = flatSeg(:,2)-1;

    % Exclude flat segments shorter than 5 samples
    smallIntervals = diff(flatSeg,[],2) <= 5;
    flatSeg(smallIntervals,:) = [];

    % Expand segments because edges that are progressively flat are missed
    if ~isempty(flatSeg)
        flatSeg(:,1) = flatSeg(:,1) - floor(EEG.srate/4);
        flatSeg(:,2) = flatSeg(:,2) + floor(EEG.srate/4);
        if flatSeg(1,1) < 0
            flatSeg(1,1) = 1;
        end
    end

    % Update flat mask for vis_artifacts
    mask = false(1,EEG.pnts);
    for iSeg = 1:size(flatSeg,1)
        lowBound = flatSeg(iSeg,1);
        highBound = flatSeg(iSeg,2);
        mask(lowBound:highBound) = true;
    end

    % Calculate total length and ratio removed (before removing them)
    flatSamples = sum(diff(flatSeg,[],2));
    flatSegRatio = round(flatSamples/EEG.pnts*100,1);
    % flatSegSec = round(flatSamples/EEG.srate,1);

    % Remove from data
    EEG = pop_select(EEG,'nopoint',flatSeg);

    % Visualize
    if vis
        EEG.etc.clean_sample_mask = ~mask;
        vis_artifacts(EEG,oriEEG);
    end

else
    flatSegRatio = 0;
    % flatSegSec = 0;
end

% Clear mask to avoid interference with ASR later
if isfield(EEG.etc,'clean_sample_mask')
    EEG.etc = rmfield(EEG.etc,'clean_sample_mask');
end
EEG = eeg_checkset(EEG);

% Print what was removed
% fprintf("%g %% were flat segments that were removed (%g min). \n", flatSegRatio, flatSegSec/60)
fprintf("%g %% were flat segments that were removed. \n", flatSegRatio)
