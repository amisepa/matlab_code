%Computes and plots means and 95% High density intervals (HDI) 
% using Bayesian bootstrap. Adapted from LIMO-EEG tools
% 
% Usage:
% plotHDI(xAxis, data1, data2, type, grp, aLevel, h, data1Name, data2Name, plot1Title)
% plotHDI(freqs, data1, data2, 'Trimmed mean', 'dependent', .05, h, 'condition1', 'condition2', 'Power spectrum (trimmed mean + 95% HDI)'); 

% data format:  e.g. frames x participants
% type: 'Mean', 'Trimmed Mean'
% grp: 'dependent', or 'independent'
% aLevel: alpha threshold probability (default = .05)
% h: true or false index for each xAxis value, to plot significance bars at the bottom
% 
% Cedric Cannard 2021


function plotHDI(xAxis, data1, data2, type, grp, alphaLevel, h, data1Name, data2Name, plot1Title)
 
fprintf('Computing Bayesian bootstrap for HDI estimation... \n')

if exist('h', 'var') || ~isempty(h)
    sigBars = true;
else
    sigBars = false;
end

%remove unecesary dimensions if any
% if length(size(data1)) > 2 && size(data1,1) == 1
%     data1 = squeeze(data1);
%     data2 = squeeze(data2);
% end

%COMPUTE HDI
if strcmp(grp, 'dependent') %%%%%%% PAIRED METHOD ONLY %%%%%%%%%%%%%
    [est1, HDI1] = computeHDI(data1, type, 1-alphaLevel);
    [est2, HDI2] = computeHDI(data2, type, 1-alphaLevel);
    [est3, HDI3] = computeHDI(data1-data2,type,1-alphaLevel);   
else
    error('Need to add independent method'); %%%%%% ADD INDPT METHOD %%%%%%%
end

% Dependent vs Independent
% if strcmp(grp, 'dependant')
%     [Ty,CI,diff,se,df,p] = yuend(est1,est2,20,alphaLevel)
% else
%     [Ty,diff,CI,da,db,df,p] = yuen(est1,est2,20,alphaLevel)
% end

%Colors
color1 = [0, 0.4470, 0.7410];           %blue
color2 = [0.8500, 0.3250, 0.0980];      %red
color3 = [0.4660, 0.6740, 0.1880];      %green


%Data1 (mean + 95% HDI)
figure; set(gcf,'Color','w');
subplot(2,1,1);
p1 = plot(xAxis, est1,'LineWidth',2,'Color', color1);
patch([xAxis fliplr(xAxis)], [HDI1(1,:) fliplr(HDI1(2,:))], ...
    color1,'FaceAlpha',.4,'EdgeColor',color1,'EdgeAlpha',0.9);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; hold on; %box on

%Data2 (mean + 95% HDI)
hold on;
p2 = plot(xAxis, est2,'LineWidth',2,'Color', color2);
patch([xAxis fliplr(xAxis)], [HDI2(1,:) fliplr(HDI2(2,:))], ...
    color2,'FaceAlpha',.4,'EdgeColor',color2,'EdgeAlpha',0.9);
% set(gca,'FontSize',12,'layer','top'); 
% grid on; axis tight; %box on
% ylabel('MSE','FontSize',12);
title([plot1Title ' (' type ' and 95% HDI)']); 

%Difference between data1 and data2 (mean + 95% HDI)
subplot(2,1,2)
plot(xAxis, est3,'LineWidth',2,'Color', color3);
patch([xAxis fliplr(xAxis)], [HDI3(1,:) fliplr(HDI3(2,:))], ...
    color3,'FaceAlpha',.4,'EdgeColor',color3,'EdgeAlpha',0.9);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; %box on
% xlabel('Frequency (Hz)','FontSize',12)
% xlabel('Time (s)','FontSize',12)
% ylabel('Difference','FontSize',12);
title(['Difference (' type ' and 95% HDI)']); 

%Plot significance bar at the bottom
if sigBars
    plotSigBar(h, xAxis);
end

legend([p1, p2], {data1Name,data2Name}, 'Orientation','vertical'); 

end    

