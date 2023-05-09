function plot_topo(Values, loc_file, figPath)

figure('color','w','Visible','off');
% fig.Visible = 'off';    % hide it immediately

% cmap = colormap('jet');
% cmaplen = size(cmap,1);
% [r,c] = size(Values);
Values = Values(:); % make Values a column vector
ContourVals = Values(:); % values for contour

% Read the channel location information
[~, labels, Th, Rd, indices] = readlocs(loc_file);
Th = pi/180*Th;    % convert degrees to radians
allchansind = 1:length(Th);
plotchans = 1:12;
plotchans = intersect_bc(plotchans, indices);

% remove infinite and NaN values
if length(Values) > 1
    inds = union_bc(find(isnan(Values)), find(isinf(Values))); % NaN and Inf values
    plotchans = setdiff_bc(plotchans, inds);
end

[x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
Th        = Th(plotchans);
Rd        = Rd(plotchans);
x         = x(plotchans);
y         = y(plotchans);
labels    = labels(plotchans);  % remove labels for electrodes without locations
labels    = strvcat(labels);    % make a label string matrix

if ~isempty(Values) && length(Values) > 1
    Values      = Values(plotchans);
    ContourVals = ContourVals(plotchans);
end

% Don't plot channels with Rd > 1 (below head)
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
% default_intrad = 1;     % indicator for (no) specified intrad
intrad = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location

% Set radius of head cartoon
rmax = 0.5;   % actual head radius
if plotrad >= rmax
    headrad = rmax;  % anatomically correct
else % if plotrad < rmax
    headrad = 0;    % don't plot head
    if strcmpi(VERBOSE, 'on')
        fprintf('cannot plot head since plotrad (%5.4g) < 0.5\n', plotrad);
    end
end

% Shrink mode
% shrinkfactor = .2;
% fprintf('Shrinking coordinates by %3.2f\n', shrinkfactor);
% plotrad = rmax/(1-shrinkfactor);
% headrad = plotrad;
% if abs(headrad-rmax) > 1e-2
%     fprintf(' Warning: With this "shrink" setting, the cartoon head will NOT be anatomically correct.\n');
% end

% Find plotting channels
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
intchans = find(x <= intrad & y <= intrad); % interpolate and plot channels inside interpolation square
% intchans = find(Rd <= intrad); % interpolate channels in the radius intrad circle only

% Eliminate channels not plotted
allx      = x;
ally      = y;
intchans; % interpolate using only the 'intchans' channels
pltchans; % plot using only indicated 'plotchans' channels
if length(Values) == length(Th)  % if as many map Values as channel locs
    intValues      = Values(intchans);
    intContourVals = ContourVals(intchans);
    Values         = Values(pltchans);
    ContourVals    = ContourVals(pltchans);
end

allchansind = allchansind(pltchans);
intTh = Th(intchans);           % eliminate channels outside the interpolation area
intRd = Rd(intchans);
intx  = x(intchans);
inty  = y(intchans);
Th    = Th(pltchans);              % eliminate channels outside the plotting area
Rd    = Rd(pltchans);
x     = x(pltchans);
y     = y(pltchans);
labels = labels(pltchans,:);

% Squeeze channel locations to <= rmax so that outermost channels will be plotted just inside rmax
squeezefac = rmax/plotrad;
intRd = intRd*squeezefac; % squeeze electrode arc_lengths towards the vertex
Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex to plot all inside the head cartoon
intx = intx*squeezefac;
inty = inty*squeezefac;
x = x*squeezefac;
y = y*squeezefac;
allx = allx*squeezefac;
ally = ally*squeezefac;

% rotate channels based on chaninfo
nosedir = '+x';
if ~strcmpi(nosedir, '+x')
    if strcmpi(nosedir, '+y')
        rotate = 3*pi/2;
    elseif strcmpi(nosedir, '-x')
        rotate = pi;
    else
        rotate = pi/2;
    end
    allcoords = (inty + intx*sqrt(-1))*exp(sqrt(-1)*rotate);
    intx = imag(allcoords);
    inty = real(allcoords);
    allcoords = (ally + allx*sqrt(-1))*exp(sqrt(-1)*rotate);
    allx = imag(allcoords);
    ally = real(allcoords);
    allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*rotate);
    x = imag(allcoords);
    y = real(allcoords);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%

% Find limits for interpolation
xmin = min(-rmax,min(intx)); xmax = max(rmax,max(intx));
ymin = min(-rmax,min(inty)); ymax = max(rmax,max(inty));

% Interpolate scalp map data
gridScale = 67;  % plot map on a 67X67 grid
xi = linspace(xmin,xmax,gridScale);   % x-axis description (row vector)
yi = linspace(ymin,ymax,gridScale);   % y-axis description (row vector)
% try
[~,~,Zi] = griddata(inty,intx,double(intValues),yi',xi,'v4'); % interpolate data
[Xi,Yi,ZiC] = griddata(inty,intx,double(intContourVals),yi',xi,'v4'); % interpolate data
% catch,
%   [Xi,Yi] = meshgrid(yi',xi);
%   Zi  = gdatav4(inty,intx,double(intValues), Xi, Yi);
%   ZiC = gdatav4(inty,intx,double(intContourVals), Xi, Yi);
% end

% Mask out data outside the head
mask = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
ii = find(mask == 0);
Zi(ii)  = NaN;                         % mask non-plotting voxels with NaNs
ZiC(ii) = NaN;                         % mask non-plotting voxels with NaNs
grid = plotrad;                       % unless 'noplot', then 3rd output arg is plotrad

% Calculate colormap limits
amax = max(max(abs(Zi)));
amin = -amax;
delta = xi(2)-xi(1); % length of grid entry

% Scale the axes
hold on
h = gca; % uses current axes
set(gca,'Xlim',[-rmax rmax]*1.3,'Ylim',[-rmax rmax]*1.3); % specify size of head axes in gca
unsh = (67+1)/67; % un-shrink the effects of 'interp' SHADING

% Set color axis
cax_sgn = sign([amin amax]);                                                  % getting sign
if ~all(cax_sgn == 0)
    %     caxis([amin+cax_sgn(1)*(0.05*abs(amin)) amax+cax_sgn(2)*(0.05*abs(amax))]);   % Adding 5% to the color limits
    clim([amin+cax_sgn(1)*(0.05*abs(amin)) amax+cax_sgn(2)*(0.05*abs(amax))]);   % Adding 5% to the color limits
end

% Plot filled ring to mask jagged grid boundary
hwidth = 0.007;                   % width of head ring
hin = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius
rwidth = 0.035;         % width of blanking outer ring
rin =  rmax*(1-rwidth/2);              % inner ring radius
if hin > rin, rin = hin;  end        % dont blank inside the head ring

% plot map and contour
tmph = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi))-0.1,Zi,...
       'EdgeColor','none','FaceColor','flat');
% contour(Xi,Yi,ZiC,6,'k');

% mask the jagged border around rmax
circGrid = 201;       % number of angles to use in drawing circles
circ = linspace(0,2*pi,circGrid);
rx = sin(circ);
ry = cos(circ);
ringx = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];
ringh = patch(ringx,ringy,0.01*ones(size(ringx)),bckgrd,'edgecolor','none'); hold on

% Plot head outline
% headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
% heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];
headx = [rx(:)' rx(1)]*hin;
heady = [ry(:)' ry(1)]*hin;
ringh = plot(headx,heady);
set(ringh,'color','k','linewidth',2); hold on

% Nose and ears parameters
base  = rmax-.0046;
basex = 0.18*rmax;                   % nose width
tip   = 1.15*rmax;
tiphw = .04*rmax;                    % nose tip half width
tipr  = .01*rmax;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = headrad/plotrad;   % squeeze the model ears and nose by this factor

% plot nose
plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
    2*ones(size([basex;tiphw;0;-tiphw;-basex])),'color','k','LineWidth',2);

% plot right ear
plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color','k','LineWidth',2)

% plot left ear
plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color','k','LineWidth',2)

% Show electrode information
plotax = gca;
axis square   
axis off
set(plotax, 'xlim', [-plotrad plotrad]);
set(plotax, 'ylim', [-plotrad plotrad]);
% pos = get(gca,'position');
% set(plotax,'position',pos);
% xlm = get(gca,'xlim');
% ylm = get(gca,'ylim');
% set(plotax,'xlim',xlm);
% set(plotax,'ylim',ylm);      
% axis equal;
% lim = [-plotrad plotrad];
% set(gca, 'xlim', [-plotrad plotrad]); 
% set(gca, 'ylim', [-plotrad plotrad]); 

% PLOT ELEC DISKS
% topoplot(dataToPlot, chanlocs, 'emarker', {'o','k',15,1});
elecHeight = 2.1;
markerSize = 15;
markerWidth = 2;
plot3(y,x,ones(size(x))*elecHeight,'o','color','k',...
    'markersize',markerSize,'linewidth',markerWidth);

% topoplot(dataToPlot, chanlocs, 'electrodes', 'labels');
% topoplot(dataToPlot,chanlocs,'electrodes', 'numpoint');

% Add values inside circles
valSize = 10;
for iChan = 1:length(Values)
    text(double(y(iChan)-0.02),double(x(iChan)),...
        elecHeight,num2str(Values(iChan)),'HorizontalAlignment','left',...
        'VerticalAlignment','middle','Color','k','FontSize',valSize, 'FontWeight','bold');
end

% colormap(blueredyellowcmap())
% colormap('parula')
colormap('jet')

colorbar; clim([0 100])
% set(gcf, 'color', bckgrd);
hold off
axis off

% save jpg image
% imwrite(fig, figPath)

%% SUBFUNCTIONS

% % Draw circle
% function h2 = disk(X, Y, radius, colorfill, oriangle, endangle, segments)
% A = linspace(oriangle/180*pi, endangle/180*pi, segments-1);
% if endangle-oriangle == 360
%  	 A  = linspace(oriangle/180*pi, endangle/180*pi, segments);
%      h2 = patch( X   + cos(A)*radius(1), Y   + sin(A)*radius(end), zeros(1,segments)+3, colorfill);
% else A  = linspace(oriangle/180*pi, endangle/180*pi, segments-1);
%      h2 = patch( [X X + cos(A)*radius(1)], [Y Y + sin(A)*radius(end)], zeros(1,segments)+3, colorfill);
% end
% set(h2, 'FaceColor', colorfill);
% set(h2, 'EdgeColor', 'none');


% % MATLAB 4 GRIDDATA interpolation
% % Reference:  David T. Sandwell, Biharmonic spline interpolation of GEOS-3 and SEASAT altimeter
% % data, Geophysical Research Letters, 2, 139-142, 1987.  Describes interpolation using value or
% % gradient of value in any dimension.
% function vq = gdatav4(x,y,v,xq,yq)
% xy = x(:) + 1i*y(:);
% d = abs(xy - xy.'); % Determine distances between points
% g = (d.^2) .* (log(d)-1);   % Determine weights for interpolation (Green's function)
% % Fixup value of Green's function along diagonal
% g(1:size(d,1)+1:end) = 0;
% weights = g \ v(:);
% [m,n] = size(xq);
% vq = zeros(size(xq));
% xy = xy.';
% % Evaluate at requested points (xq,yq).  Loop to save memory.
% for i=1:m
%     for j=1:n
%         d = abs(xq(i,j) + 1i*yq(i,j) - xy);
%         g = (d.^2) .* (log(d)-1);   % Green's function.
%
%         % Value of Green's function at zero
%         g(d==0) = 0;
%         vq(i,j) = g * weights;
%     end
% end


