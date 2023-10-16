%Computes and plots means and 95% High density intervals (HDI) 
% using Bayesian bootstrap. Adapted from LIMO-EEG tools
% 
% Usage:
% plotHDI(xAxis, data1, data2, type, grp, alphaLevel, h, data1Name, data2Name, plot1Title)
% plotHDI(freqs, data1, data2, 'Trimmed mean', 'dependent', .05, h, 'condition1', 'condition2', 'Power spectrum (trimmed mean + 95% HDI)'); 

% data format:  e.g. frames x participants
% type: 'Mean', 'Trimmed Mean'
% grp: 'dependent', or 'independent'
% aLevel: alpha threshold probability (default = .05)
% h: true or false index for each xAxis value, to plot significance bars at the bottom
% 
% Cedric Cannard 2021


function plotHDI(xAxis, data1, data2, type, grp, alphaLevel, h, data1Name, data2Name)
 
fprintf('Computing Bayesian bootstrap for HDI estimation... \n')

%remove unnecessary dimensions if any
% if length(size(data1)) > 2 && size(data1,1) == 1
%     data1 = squeeze(data1);
%     data2 = squeeze(data2);
% end

% COMPUTE 95% high-density intervals (HDI)
if strcmp(grp, 'dependent') %%%%%%% PAIRED METHOD ONLY %%%%%%%%%%%%%
    [est1, HDI1] = computeHDI(data1, type, 1-alphaLevel);
    [est2, HDI2] = computeHDI(data2, type, 1-alphaLevel);
    [est3, HDI3] = computeHDI(data1-data2,type,1-alphaLevel);   
else
    error('Need to add independent method'); %%%%%% ADD INDPT METHOD %%%%%%%
end

% Colors
color1 = [0, 0.4470, 0.7410];           % blue
color2 = [0.8500, 0.3250, 0.0980];      % red
color3 = [0.4660, 0.6740, 0.1880];      % green
% cb = cbrewer2('qual', 'Set3', 12, 'pchip');
% color1 = cb(5,:);  % 5=blue, 4=red, 1=green, 2=yellow
% color2 = cb(4,:);  % red
% color3 = cb(1,:);  % green

% Data1 (mean + 95% HDI)
figure('color','w'); 
subplot(2,1,1);
p1 = plot(xAxis,est1,'LineWidth',1,'Color', color1);
patch([xAxis fliplr(xAxis)], [HDI1(1,:) fliplr(HDI1(2,:))], ...
    color1,'FaceAlpha',.3,'EdgeColor',color1,'EdgeAlpha',0.9);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; hold on; box on

% Data2 (mean + 95% HDI)
hold on;
p2 = plot(xAxis, est2,'LineWidth',1,'Color', color2);
patch([xAxis fliplr(xAxis)], [HDI2(1,:) fliplr(HDI2(2,:))], ...
    color2,'FaceAlpha',.3,'EdgeColor',color2,'EdgeAlpha',0.9);
title(sprintf('%s + %g%% HDI',type,(1-alphaLevel)*100)); 

% Difference between data1 and data2 (mean + 95% HDI)
subplot(2,1,2)
plot(xAxis, est3,'LineWidth',1,'Color', color3);
patch([xAxis fliplr(xAxis)], [HDI3(1,:) fliplr(HDI3(2,:))], ...
    color3,'FaceAlpha',.3,'EdgeColor',color3,'EdgeAlpha',0.9);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; box on
% title(sprintf('Difference (%s + 95% HDI)',type)); 
ylabel('Difference','FontSize',11,'FontWeight','bold')

% Plot significance bar at the bottom
if isempty(h)
    for i = 1:length(xAxis)
        if HDI3(1,i)<0 && HDI3(2,i)<0 || HDI3(1,i)>0 && HDI3(2,i)>0
            h(i) = true;
        else
            h(i) = false;
        end
    end
end
plotSigBar(h, xAxis);

legend([p1, p2], {data1Name,data2Name}, 'Orientation','vertical'); 
set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

