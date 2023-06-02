%% Get the number of principal components that should be retained to
% explain X% of variance in the data. This is useful for
% diemnsion-reduction to use the minimum number of components to retain. 
% 
% Inputs:
%   data
%   varThresh   - variance threshold to retain in % (default = 99)
%   vis         - visualize (true) or not (false)
% 
% Usage: 
%   nComps = get_nPCA(data, varThresh, vis)
%   nComps = get_nPCA(data, 99, true)
% 
% Copyright (C) - Cedric Cannard, 2023

function nComps = get_nPCA(data, varThresh, vis)

if ~exist('varThresh','var') || isempty(varThresh)
    varThresh = 99;
end
if ~exist('varThresh','var') || isempty(vis)
    vis = 1;
end

% Calculate cumulative variance explained
[coeff, score, ~, ~, explained] = pca(double(data)); 
cumulativeExplained = cumsum(explained);

% number of components explaining 99% of variance
nComps = find(cumulativeExplained >= varThresh, 1); 

% Plot cumulative variance explained
if vis
    figure('color','w');
    plot(cumulativeExplained, 'o-');
    xlabel('Number of Components');
    ylabel('Cumulative Variance Explained (%)');
    title(sprintf('# of components to retain to preserve 99%% of variance = %g', nComps));
    grid on; ylim([0 100])
end
