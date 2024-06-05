function EEG = clean_eeg(EEG,thresh,baseline,vis)

if isfield(EEG.etc,'clean_sample_mask')
    EEG.etc = rmfield(EEG.etc,'clean_sample_mask');
end

oriEEG2 = EEG;      % copy for vis_artifacts

% Identify bad segments with ASR
m = memory; maxmem = round(.5*(m.MemAvailableAllArrays/1000000),1); % use 50% of available RAM
cleanEEG = clean_asr(EEG,thresh,[],[],[],baseline,[],[],false,false,maxmem);

% Create mask
mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
EEG.etc.clean_sample_mask = true(1, length(mask)); % initialize all samples as clean
badData = reshape(find(diff([false mask false])), 2, [])';
badData(:, 2) = badData(:, 2)-1;


if ~isempty(badData)

    % exclude very short segments < 10 samples
    smallIntervals = diff(badData')' < 10;  
    badData(smallIntervals, :) = [];
    for i = 1:size(badData, 1)
        EEG.etc.clean_sample_mask(badData(i, 1):badData(i, 2)) = false;
    end
    
    % Remove them
    EEG = pop_select(EEG,'nopoint',badData);
    fprintf('%g %% of data were considered to be artifacts and were removed. \n', round((1-EEG.xmax/oriEEG2.xmax)*100,2))

    % Visualize what was removed
    if vis
        vis_artifacts(EEG,oriEEG2,'ShowSetname',false);
        set(gcf,'Name','Artifacts removed by ASR','NumberTitle','Off','Toolbar','none','Menu','none','color','w'); box on;
    end
end
