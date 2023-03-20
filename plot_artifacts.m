function plot_artifacts(cleanEEG, rawEEG)


sample_mask = sum(abs(rawEEG.data-cleanEEG.data),1) < 1e-8;
retain_data_intervals = reshape(find(diff([false sample_mask false])),2,[])';
retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;
% if ~isempty(retain_data_intervals)
%     smallIntervals = diff(retain_data_intervals')' < 10;
%     retain_data_intervals(smallIntervals,:) = [];
% end
% fprintf('%3.1f %% of data were removed. \n', cleanEEG.xmax/rawEEG.xmax*100)

cleanEEG = pop_select(cleanEEG, 'point', retain_data_intervals);
cleanEEG.etc.clean_sample_mask = sample_mask;
vis_artifacts(cleanEEG,rawEEG);


% Mask
% mask = sum(abs(rawEEG.data-cleanEEG.data),1) > 1e-10;
% rawEEG.etc.clean_sample_mask = ~mask;
% clear cleanEEG;
% 
% % Bad data segments
% badData = reshape(find(diff([false mask false])),2,[])';
% badData(:,2) = badData(:,2)-1;
% 
% % remove segments shorter than 10 samples
% if ~isempty(badData)
%     smallIntervals = diff(badData')' < 10;
%     badData(smallIntervals,:) = [];
% end
% 
% % Flatten them for comparison
% cleanEEG = rawEEG;
% for iSegment = 1:size(badData,1)
%     cleanEEG.data(:,badData(iSegment,1):badData(iSegment,2)) = 0;
% end

% EEG = pop_select(EEG,'nopoint',badData);

% vis_artifacts(cleanEEG,rawEEG);  % to compare (doesn't work)

