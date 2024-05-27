% Assuming data1 and data2 are your two conditions with size [channels, time points, subjects]
% data1 and data2 are matrices of size [channels, time points, subjects]

% data1 = randn(1, 750, 78); % Example data, replace with your actual data
% data2 = randn(1, 750, 78); % Example data, replace with your actual data
% 
% % parameters
% nPerm = 1000;
% alphalvl = 0.05;

function [tvals, pvals] = run_perm_stats(data1,data2,alphalvl,nPerm)

% Add an extra dimension if needed (one channel data squeezed)
if ndims(data1) == 2 
    data1 = reshape(data1, 1, size(data1, 1), size(data1, 2));
    data2 = reshape(data2, 1, size(data2, 1), size(data2, 2));
end

% Initialize variables
nChan = size(data1, 1);
nTimes = size(data1, 2);
nSub = size(data1, 3);
tvals = zeros(nChan, nTimes);           % Observed t-values
permDist = zeros(nChan, nTimes, nPerm); % Permutation distribution
pvals = zeros(nChan, nTimes);           % p-values from permutation test

% Observed t-values
disp("Running statistical tests on observed data (all electrodes)")
for iChan = 1:nChan
    fprintf('')
    parfor t = 1:nTimes
        if strcmpi(dpt, 'dpt')
            [tval,~,~,~,~,pval] = yuend(squeeze(data1(iChan,t,:)),squeeze(data2(iChan,t,:)),20,alphalvl);   % paired for 2D vector
        else
            [tval,~,~,~,~,pval] = yuen(squeeze(data1(iChan,t,:)),squeeze(data2(iChan,t,:)),20,alphalvl);   % unpaired for 2D vector
        end
        % [~, ~, ~, stats] = ttest(squeeze(data1(iChan, t, :)), squeeze(data2(iChan, t, :)));
        tvals(iChan, t) = tval;
    end
end

% Permutation test
for perm = 1:nPerm
    permData1 = data1;
    permData2 = data2;
    
    parfor subj = 1:nSub
        if rand > 0.5
            permData1(:, :, subj) = data2(:, :, subj);
            permData2(:, :, subj) = data1(:, :, subj);
        end
    end
    
    for iChan = 1:nChan
        parfor t = 1:nTimes
            [~, ~, ~, permStats] = ttest(squeeze(permData1(iChan, t, :)), squeeze(permData2(iChan, t, :)));
            permDist(iChan, t, perm) = permStats.tstat;
        end
    end
end

% Calculate p-values
for iChan = 1:nChan
    parfor t = 1:nTimes
        pvals(iChan, t) = mean(abs(permDist(iChan, t, :)) >= abs(tvals(iChan, t)));
    end
end

