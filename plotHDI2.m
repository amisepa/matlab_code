% Computes and plots trimmed means and 95% High density intervals (HDI) 
% using Bayesian bootstrap. Adapted from LIMO-EEG tools. 
% 
% plotHDI2 only does upper plot, no difference plot below.
% 
% Usage:
% plotHDI2(xAxis, data1, data2, type, grp, alphaLevel, h, data1Name, data2Name, plot1Title)% 
% type: 'Mean', 'Trimmed Mean'
% grp: 'dependent', or 'independent'
% alphaLevel: predetermined threshold probability (e.g. 0.05)
% h: true or false index for each xAxis value, to plot significance bars at the bottom
% 
% Cedric Cannard 2021


function plotHDI2(xAxis, data1, data2, type, grp, alphaLevel, h, data1Name, data2Name, plot1Title)
 
fprintf('Computing Bayesian bootstrap for HDI estimation... \n')

if exist('h', 'var') || ~isempty(h)
    sigBars = true;
else
    sigBars = false;
end

%remove unnecessary dimensions if any
% if length(size(data1)) > 2 && size(data1,1) == 1
%     data1 = squeeze(data1);
%     data2 = squeeze(data2);
% end

%COMPUTE HDI
if strcmp(grp, 'dependent') %%%%%%% PAIRED METHOD ONLY %%%%%%%%%%%%%
    [est1, HDI1] = computeHDI(data1, type, 1-alphaLevel);
    [est2, HDI2] = computeHDI(data2, type, 1-alphaLevel);
else
    error('Need to add independent method'); %%%%%% ADD INDPT METHOD %%%%%%%
end

%Colors
color1 = [0, 0.4470, 0.7410];           %blue
color2 = [0.8500, 0.3250, 0.0980];      %red

%PLOT MEAN AND 95% HDI OF EACH CONDITION/GROUP
% figure; set(gcf,'Color','w');
p1 = plot(xAxis, est1,'LineWidth',2,'Color',color1); hold on;
fillhandle = fill([xAxis fliplr(xAxis)],[HDI1(1,:) fliplr(HDI1(2,:))], color1);
set(fillhandle,'EdgeColor', color1,'FaceAlpha',0.3,'EdgeAlpha',0.8);
p2 = plot(xAxis, est2,'LineWidth',2,'Color', color2);
fillhandle = fill([xAxis fliplr(xAxis)],[HDI2(1,:) fliplr(HDI2(2,:))], color2);
set(fillhandle,'EdgeColor', color2,'FaceAlpha',0.3,'EdgeAlpha',0.8);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; box on
% ylabel('PSD (dB)','FontSize',12);
% xlabel('Frequency (Hz)','FontSize',12);
% ylabel('Alpha asymmetry amplitude','FontSize',12);
% title([plot1Title ' (' type ' and 95% HDI)']); 
title(plot1Title); 
% hPlots = flip(findall(gcf,'Type','Line')); % flipped, because the lines are found in reverse order of appearance.
% legend(hPlots, {data1Name,data2Name}, 'Orientation','vertical'); 

%Plot significance bar at the bottom
if sigBars
    plotSigBar(h, xAxis);
end

legend([p1, p2], {data1Name,data2Name}, 'Orientation','vertical'); 
    

