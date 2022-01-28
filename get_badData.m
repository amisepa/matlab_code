%% Outputs segments of bad data removed by the clean_rawdata plugin
% 
% Example: badData = getbadData(EEG)
% 
% Cedric Cannard, 2020

function bad_data = get_badData(EEG)

mask = ~EEG.etc.clean_sample_mask;
bad_data = [];
count = 1;
for a = 1:length(mask)
    
    %For last sample: end of bad period
    if mask(a) == 1 && a == length(mask)
        bad_data(count,2) = a;
        break
    end
    
    %For 1st sample: beginning of bad period
    if mask(a) == 1 && a == 1
        bad_data(count,1) = a;
        continue
    end
    
    %For other samples: beginning of bad period
    if mask(a) == 1 && mask(a-1) == 0 
        bad_data(count,1) = a;
    end
    
    %For other samples: end of bad period
    if mask(a) == 1 && mask(a+1) == 0
        bad_data(count,2) = a;
        count = count+1;
    end
end

%if space between 2 bad epochs is very small, merge them into one epoch
spacer = 32;
for a = 1:size(bad_data,1)-1
    if bad_data(a,2) + spacer >= bad_data(a+1,1)
        bad_data(a,2) = bad_data(a+1,1);
    end
end
bad_data2 = bad_data;
for a = 1:size(bad_data2,1)
    if a == size(bad_data2,1)
        break
    end
    section = bad_data2(a,1) : bad_data2(a,2);
    section_next = bad_data2(a+1,1) : bad_data2(a+1,2);
    if sum(ismember(section, section_next))
        bad_data2(a,1) = min(bad_data2(a,1), bad_data2(a+1,1));
        bad_data2(a+1,1) = min(bad_data2(a,1), bad_data2(a+1,1));
        bad_data2(a,2) = max(bad_data2(a,2), bad_data2(a+1,2));
        bad_data2(a+1,2) = max(bad_data2(a,2), bad_data2(a+1,2));
    end
end
bad_data2 = unique(bad_data2, 'rows');    
bad_data = bad_data2;
bad_data2 = [];

%Do this a second time
spacer = 32;
for a = 1:size(bad_data,1)-1
    if bad_data(a,2) + spacer >= bad_data(a+1,1)
        bad_data(a,2) = bad_data(a+1,1);
    end
end
bad_data2 = bad_data;
for a = 1:size(bad_data2,1)
    if a == size(bad_data2,1)
        break
    end
    section = bad_data2(a,1) : bad_data2(a,2);
    section_next = bad_data2(a+1,1) : bad_data2(a+1,2);
    if sum(ismember(section, section_next))
        bad_data2(a,1) = min(bad_data2(a,1), bad_data2(a+1,1));
        bad_data2(a+1,1) = min(bad_data2(a,1), bad_data2(a+1,1));
        bad_data2(a,2) = max(bad_data2(a,2), bad_data2(a+1,2));
        bad_data2(a+1,2) = max(bad_data2(a,2), bad_data2(a+1,2));
    end
end
bad_data2 = unique(bad_data2, 'rows');    
bad_data = bad_data2;
bad_data2 = [];

%get duration of each epoch of bad data
for b = 1:size(bad_data,1)
    bad_data(b,3) = bad_data(b,2)+1 - bad_data(b,1);
end

%% Old method to get bad data after manual cleaning
% if isempty(EEG.event)   %if no bad samples at all, remove the first one
%     badData = [1 2];
% else
%     for n = 1:length(EEG.event)
%         badData(n,:) = [round(EEG.event(n).latency) round(EEG.event(n).latency)-1+EEG.event(n).duration];
%         if n > 1
%             tmp = badData(n-1,2) + EEG.event(n).latency - EEG.event(n-1).latency;
%             badData(n,:) = [tmp tmp + EEG.event(n).duration];
%         end
%     end
%     badData = sort(badData);
% end
