%% Detect flat lines in EEG data, excludes very small segments (< 5 samples),
% plots them in red (superimposed on raw EEG data), and print portion of data
% that were removed (seconds and ratio).
%
% Copyright (C) - Cedric Cannard, 2023

function EEG = rm_flatSegments(EEG,flatThresh,vis)

% Threshold
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
flat_mask = nan(nChan,EEG.pnts+1);
for iChan = 1:nChan
    flat_mask(iChan,:) = [false abs(diff(EEG.data(iChan,:)))<flatThresh false]; % false 1st and last sample to preserve original size
    % flat_mask(iChan,:) = [false isoutlier(abs(diff(EEG.data(iChan,:))),'mean') false];   % this can detect EEG artifacts
end
if nChan == 4
    flat_mask_front = or(flat_mask(2,:),flat_mask(3,:));
    flat_mask_post = or(flat_mask(1,:),flat_mask(4,:));
    flat_mask = or(flat_mask_front, flat_mask_post);
elseif nChan == 3
    flat_mask_tmp = or(flat_mask(1,:),flat_mask(2,:));
    flat_mask = or(flat_mask(3,:),flat_mask_tmp);
elseif nChan == 2
    flat_mask = or(flat_mask(1,:),flat_mask(2,:));
end
flat_mask(end) = [];

% Remove them if present
if sum(flat_mask) > 0
    flatSeg = reshape(find(diff([false flat_mask false])),2,[])';
    flatSeg(:,2) = flatSeg(:,2)-1;

    % Exclude flat segments shorter than 5 samples
    smallIntervals = diff(flatSeg,[],2) <= 5;
    % for iSeg = 1:size(flatSeg,1) % update mask
    %     if smallIntervals(iSeg)
    %         lowBound = flatSeg(iSeg,1);
    %         highBound = flatSeg(iSeg,2);
    %         flat_mask(lowBound:highBound) = false;
    %     end
    % end
    flatSeg(smallIntervals,:) = [];

    % Expand segments because edges that are progressively flat are missed
    flatSeg(:,1) = flatSeg(:,1) - floor(EEG.srate/4);
    flatSeg(:,2) = flatSeg(:,2) + floor(EEG.srate/4);
    if flatSeg(1,1) < 0
        flatSeg(1,1) = 1;
    end

    % Update flat mask for vis_artifacts
    flat_mask = false(1,EEG.pnts);
    for iSeg = 1:size(flatSeg,1)
        lowBound = flatSeg(iSeg,1);
        highBound = flatSeg(iSeg,2);
        flat_mask(lowBound:highBound) = true;
    end

    % Remove from data
    EEG = pop_select(EEG,'nopoint',flatSeg);

    % Visualize
    if vis
        EEG.etc.clean_sample_mask = ~flat_mask;
        vis_artifacts(EEG,oriEEG);
    end

    % EEG = cleanEEG2;

    % calculate total length and ratio removed
    flatSamples = sum(diff(flatSeg,[],2));
    flatSegRatio = round(flatSamples/EEG.pnts,2);
    flatSegSec = round(flatSamples/EEG.srate,1);

else
    flatSegRatio = 0;
    flatSegSec = 0;
end

% Print what was removed
fprintf("%g %% were flat segments that were removed (%g min). \n", flatSegRatio, flatSegSec/60)
