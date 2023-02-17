% Get latency bounds of artifactual segments and merge them if gap
% separating them is very small to reduce number of discontinuities.
% 
% Inputs:
%   - clean EEG structure
%   - raw EEG structure
% 
% Output:
%   - bad data segments
% 
% Cedric Cannard, 2022

function [bad_data, mask] = get_badSegments(newEEG,oldEEG)

% Portion of data which have changed
mask = sum(abs(oldEEG.data-newEEG.data),1) < 1e-10;
mask = ~mask;

% Latency bounds of bad segments
bad_data = reshape(find(diff([false mask false])),2,[])';
bad_data(:,2) = bad_data(:,2)-1;
