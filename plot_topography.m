% Plots scalp topography of MiniQ metrics using spherical spline
% interpolation. 
% 
% Inputs:
%   dataToPlot  - metric to plot (1 x nChan vector)
%   chanlocs    - structure with channel labels and coordinates

function plot_topography(dataToPlot,chanlocs)

% Parameters
shrink = .7;
radius = .5;
headWidth = 2;
contour = 'off';
elecSize = 22;
elecWidth = 2;
txtSize = 12;

% Open figure in background
% fig = figure('color','white','visible','off');

% % Score and data type
% if ~contains(figPath, 'total')
%     dataType = extractBetween(figPath,'_','_');
%     dataType = dataType{:};
% else
%     dataType = extractBetween(figPath,'_','.');
%     dataType = dataType{:};
% end
% score = extractBetween(figPath, fileparts(figPath), '_');
% score = score{:}(2:end);

% % Colormap
% if strcmpi(cmap, 'miniQmap')
%     cmap = load('-ascii', fullfile('ressources', 'miniqcolormap.txt'));
% elseif strcmpi(cmap, 'yellowredblue')
%     cmap = yellowredbluecmap;
% elseif strcmpi(cmap, 'coolhot')
%     cmap = coolhotcmap;
% elseif strcmpi(cmap, 'blueredyellow')
%     cmap = blueredyellowcmap;
% elseif strcmpi(cmap, 'bluered')
%     cmap = redbluecmap;
%     cmap = cmap(end:-1:1,:);
% elseif strcmpi(cmap, 'redblue')
%     cmap = redbluecmap;
% elseif strcmpi(cmap, 'jet2')
%     cmap = jet(128); 
%     cmap(1:64,:) = [];
% elseif strcmpi(cmap, 'hsv2')
%     cmap = hsv(128); 
%     cmap(1:64,:) = [];
% else
%     fprintf('colormap name not recognized, using default miniQmap. \n');
%     cmap = load('-ascii', fullfile('ressources', 'miniqcolormap.txt'));
% end    

pnts = linspace(0,2*pi,200);
xx = sin(pnts)*radius;
yy = cos(pnts)*radius;

% spherical plotting
xelec = [chanlocs.X];
yelec = [chanlocs.Y];
zelec = [chanlocs.Z];
dist = sqrt(xelec.^2+yelec.^2+zelec.^2);
xelec = xelec./dist;
yelec = yelec./dist;
zelec = zelec./dist;

% Shrink
[th, phi, rad] = cart2sph(xelec, yelec, zelec);
phi = (phi-pi/2)*shrink+pi/2;
[xelec, yelec, zelec] = sph2cart(th, phi, rad);

% Spherical spline interpolation
sphRes = 20;
[xsph,ysph,zsph] = sphere(sphRes);
xsph(1:(length(xsph)-1)/2,:) = [];
ysph(1:(length(xsph)-1)/2,:) = [];
zsph(1:(length(xsph)-1)/2,:) = [];
Gelec = compute_g(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = compute_g(xsph,ysph,zsph,xelec,yelec,zelec);
meanvalues = mean(dataToPlot);
values = dataToPlot - meanvalues; % zero-mean
C = pinv([Gelec;ones(1,length(Gelec))]) * [values(:);0];
valsph = zeros(1,size(Gsph,1));
for j = 1:size(Gsph,1)
    valsph(j) = sum(C .* Gsph(j,:)');
end
valsph = valsph + meanvalues;
valsph = reshape(valsph, size(xsph));

% Plot colormap of scores
topVal = max(abs(valsph(:)))*1000;
surf(-ysph/2,xsph/2,zsph/2, double(valsph), 'edgecolor', 'none'); view([0 0 1]);hold on;
shading interp;

if strcmpi(contour, 'on')
    contour3(-ysph/2, xsph/2, valsph, 5); view([0 0 1]);
end

% Coordinates for electrodes
xelec(zelec < 0) = [];
yelec(zelec < 0) = [];
x = yelec/2;
y = xelec/2;

% Plot head circle
plot3(xx,yy,ones(size(xx))*topVal, 'k', 'linewidth', headWidth); hold on;

% Plot ears & nose
earx  = [0.4960    0.505    0.520    0.530    0.540    0.533   0.550    0.543    0.530   0.500    0.490]; % rmax = 0.5
eary  = [0.0655    0.0855   0.086    0.082    0.066    0.015   -0.073   -0.09   -0.11   -0.115   -0.1];
plot3(earx,eary,ones(size(earx))*topVal,'color','k','LineWidth',headWidth)    % plot left ear
plot3(-earx,eary,ones(size(earx))*topVal,'color','k','LineWidth',headWidth)   % plot right ear
nosex = [0.0900 0.0400 0.010 0 -0.010 -0.0400 -0.0900];
nosey = [0.492  0.5300 0.5700 0.570 0.5700 0.5300 0.492 ];
plot3(nosex, nosey, ones(size(nosex))*topVal,'color','k','LineWidth',headWidth)   % plot right ear

% Plot electrodes circles
elecHeight(1:length(x)) = 2.1; %2.1 double(topVal)
% rad = sqrt(x.^2 + y.^2);
% x(rad > 0.5) = []; 
% y(rad > 0.5) = [];
plot3(-x,y,elecHeight,'o','color','k','markersize',elecSize,'linewidth',elecWidth);

% Fill color inside circles
% mapval = round((dataToPlot*64)/100);
% mapval(mapval <= 0) = 1;
% diskcolor = cmap(mapval,:);
% for iChan = 1:size(dataToPlot,2)
%     plot3(-x(iChan),y(iChan), elecHeight(iChan),'o','color','k',...
%         'MarkerFaceColor',diskcolor(iChan,:),'markersize',elecSize,...
%         'linewidth',elecWidth);
% end
diskcolor = 'w';
plot3(-x,y,elecHeight,'o','color','k','MarkerFaceColor',diskcolor, ...
    'markersize',elecSize,'linewidth',elecWidth);

% Offset and rounding
% switch score
%     case 'raw'
%         switch dataType
%             case 'abs'
%                 dataToPlot = round(dataToPlot);
%                 offset = -0.0190;
%             case 'rel'
%                 dataToPlot = round(dataToPlot,2);
%                 offset = -0.06;
%         end
% 
%     case 'percentiles'
        % offset = -0.04;
% end
offset = 0.02;

% Add values inside circles
for iChan = 1:length(dataToPlot)
%     text(-x(iChan)+offset, y(iChan), elecHeight(iChan), num2str(dataToPlot(iChan)), ...
%       'color','k','fontweight','bold','fontsize',txtSize);
    % if dataToPlot(iChan) <= 20 || dataToPlot(iChan) >= 80
        % tcolor = 'w';
    % else
        % tcolor = 'k';
    % end

    % adjust offset depending on score type

    % if strcmpi(score,'raw')
        if dataToPlot(iChan) < 100
            h = text(-x(iChan)+offset/2.1,y(iChan),elecHeight(iChan),num2str(dataToPlot(iChan)));
        else
            h = text(-x(iChan)+(offset/2.1-0.01),y(iChan),elecHeight(iChan),num2str(dataToPlot(iChan)));
        end
    % else 
    %     if dataToPlot(iChan) < 100
    %         h = text(-x(iChan)+offset,y(iChan),elecHeight(iChan),num2str(dataToPlot(iChan)));
    %     else
    %         h = text(-x(iChan)+(offset-0.01),y(iChan),elecHeight(iChan),num2str(dataToPlot(iChan)));
    %     end
    % end
    set(h,'color','k','fontweight','bold','fontsize',txtSize);
end

% Colormap and limits
% colormap(cmap);
% switch score
%     case 'percentile'
%         clim([0 100]);
%     case 'zscore'
%         clim([-2 2]);
%     case 'raw'
%         try 
clim([min(dataToPlot) max(dataToPlot)])
%         catch
%             pause
%         end
% end

% clean up
% axis([-1 1 -0.6 0.6]);
% axis tight
axis off; 
% set(gca,'ydir','normal');

% Save and close
% print(fig,'-djpeg',sprintf('%s',figPath),'-r150'); % 150 dpi
% close(fig);                   

%% Compute G for spherical spline interpolation

function g = compute_g(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - ((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +...
    (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
    (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2)/2;
g = zeros(length(x(:)),length(xelec));
m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    L = legendre(n,EI);
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);
