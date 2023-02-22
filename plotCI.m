%% Plot histogram of 2 variables and their 95% CIs
% 
% Cedric Cannard, Feb 2023

function plotCI(x1,x2,binsize)


figure('color','w');
goodCI(1) = prctile(x1, (100-95)/2);
goodCI(2) = prctile(x1, 100-(100-95)/2);
[y,x] = histcounts(x1,binsize);
y = y./max(y); x = (x(1:end-1)+x(2:end))/2;
bar(x,y,'facecolor','blue','EdgeColor','b','facealpha',.8); hold on
patch(goodCI([1 1 2 2]), [0 1.5 1.5 0],'blue', 'facealpha',.3,'edgecolor','none')
plot([1 1]*mean(x1), [0 1.5],'b--', 'linew',2, 'linewidth', 3)
badCI(1) = prctile(x2, (100-95)/2);
badCI(2) = prctile(x2, 100-(100-95)/2);
[y,x] = histcounts(x2,binsize);
y = y./max(y); x = (x(1:end-1)+x(2:end))/2;
bar(x,y,'facecolor','r','EdgeColor','r','facealpha',.8); hold on
patch(badCI([1 1 2 2]), [0 1.5 1.5 0],'red', 'facealpha',.3,'edgecolor','none')
plot([1 1]*mean(x2), [0 1.5],'r--', 'linew',2, 'linewidth', 3)
legend('x1', 'x1 (95% CI)', 'x1 (mean)', 'x2 (data)', 'x2 (95% CI)', 'x2 (mean)')
