function EEG = apply_car(EEG)
% APPLY_CAR - Applies a modified common average reference (CAR) while preserving data rank.
%
% This function adds a zero-filled surrogate for the initial reference channel, applies 
% average re-referencing that preserves effective rank, and removes the surrogate afterward.
%
% Motivation:
% Standard average referencing in EEGLAB (prior to v2023) can cause *effective rank deficiency* 
% if the initial reference is not included. This may produce ghost components in ICA.
% This implementation follows the recommendations from:
%
%    Kim, H., Luo, J., Chu, S., Cannard, C., Hoffmann, S., & Miyakoshi, M. (2023).
%    ICA's bug: How ghost ICs emerge from effective rank deficiency caused by EEG 
%    electrode interpolation and incorrect re-referencing. *Frontiers in Signal Processing*, 3:1064138.
%    https://doi.org/10.3389/frsip.2023.1064138
%
% The method:
%   1. Appends a zero-filled dummy channel labeled 'initialReference'.
%   2. Applies average referencing using EEGLAB's pop_reref(), with the dummy included.
%   3. Removes the dummy to restore the original channel count.
%
% Input:
%   EEG - EEGLAB EEG structure
%
% Output:
%   EEG - EEG structure re-referenced to average using modified CAR that preserves rank

% 1) Add a zero-filled surrogate for the initial reference
refLabel = 'initialReference';
EEG.data(end+1, :) = 0;
EEG.nbchan = EEG.nbchan + 1;
EEG.chanlocs(end+1).labels = refLabel;  % minimal fields sufficient for pop_reref

% 2) Apply average reference, explicitly setting refloc for traceability and metadata
EEG = pop_reref(EEG, [], 'refloc', struct('labels', refLabel));

% 3) Remove the dummy channel to return to the original channel count
EEG = pop_select(EEG, 'nochannel', {refLabel});
end
