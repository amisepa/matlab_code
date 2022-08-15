% Plots 2 times series, their 95% CI and significance bars at the bottom
% from h vector (FDR-corrected p-values).
% 
% Usage:
%       plotDiff(xAxis, data1, data2, h, data1Name, data2Name)
% 
% Data must be 2-D. Values in column 1 and subjects in column 2 (e.g.,freqs x subjects)
% 
% Cedric Cannard, 2021

function plotDiff(xAxis, data1, data2, h, data1Name, data2Name)

if size(xAxis,1)>size(xAxis,2)
    xAxis = xAxis';
end

if length(h) ~= length(xAxis)
    error('x-axis and h have different length.');
end
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];

nSubj = size(data1,2);
data1_mean = mean(data1,2,'omitnan');
% data1_mean = trimmean(data1,20,2);
data1_se = std(data1,[],2,'omitnan') ./ sqrt(nSubj)';      %Standard error
data1_t = tinv([.025 .975],nSubj-1);  %t-score
data1_CI = data1_mean' + (-data1_t.*data1_se)';

nSubj = size(data2,2);
data2_mean = mean(data2,2,'omitnan');
% data2_mean = trimmean(data2,20,2);
data2_se = std(data2,[],2,'omitnan') ./ sqrt(nSubj)';      %Standard error
data2_t = tinv([.025 .975],nSubj-1);  %t-score
data2_CI = data2_mean' + (-data2_t.*data2_se)';

% figure; set(gcf,'Color','w');

%Data1 (mean + 95% CI)
p1 = plot(xAxis,data1_mean,'LineWidth',2,'Color', color1);
patch([xAxis fliplr(xAxis)], [data1_CI(1,:) fliplr(data1_CI(2,:))], ...
    color1,'FaceAlpha',.4,'EdgeColor',color1,'EdgeAlpha',.7);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; hold on; %box on

%Data2 (mean + 95% CI)
p2 = plot(xAxis,data2_mean,'LineWidth',2,'Color', color2);
patch([xAxis fliplr(xAxis)], [data2_CI(1,:) fliplr(data2_CI(2,:))], ...
    color2,'FaceAlpha',.4,'EdgeColor',color2,'EdgeAlpha',.7);
set(gca,'FontSize',12,'layer','top'); 
% hold off;

%Plot significance bar at the bottom
if exist('h', 'var') || ~isempty(h)
    plotSigBar(h, xAxis);
end

legend([p1, p2], {data1Name,data2Name}, 'Orientation','vertical'); 

end
