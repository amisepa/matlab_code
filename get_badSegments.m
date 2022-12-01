% Get latency bounds of artifactual segments
% 
% Inputs:
%   - mask: array of zeros (good sample) and ones (bad samples)
%   - art_zise: minimum size of artifactiual segment to be considered artifact
%   - 

function out = get_badSegments(mask, art_size)

out = reshape(find(diff([false mask false])),2,[])';
out(:,2) = out(:,2)-1;
dur = out(:,2) - out(:,1);
out(dur<art_size,:) = [];
