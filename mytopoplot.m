function mytopoplot(chanlocs, data)

chanlocs = EEG.chanlocs;
rmax = 0.5;
rotate = 0;
EMARKERSIZE = 10;
MAPLIMITS = 'absmax';   % absmax, maxmin, [values]
ELECTRODES = 'on';                     % do
STYLE = 'blank';

% plotgrid = 'off';
% plotchans = [];
% noplot  = 'off';
% handle = [];
% Zi = [];
% chanval = NaN;
rmax = 0.5;             % actual head radius - Don't change this!
INTERPLIMITS = 'head';  % head, electrodes
INTSQUARE = 'on';       % interpolate electrodes located though the whole square containing the plotting disk
default_intrad = 1;     % indicator for (no) specified intrad
GRID_SCALE = 67;        % plot map on a 67X67 grid
CIRCGRID   = 201;       % number of angles to use in drawing circles
AXHEADFAC = 1.3;        % head to axes scaling factor
CONTOURNUM = 6;         % number of contour levels to plot
STYLE = 'both';         % default 'style': both,straight,fill,contour,blank
HEADCOLOR = [0 0 0];    % default head color (black)
CCOLOR = [0.2 0.2 0.2]; % default contour color
ELECTRODES = [];        % default 'electrodes': on|off|label - set below
MAXDEFAULTSHOWLOCS = 64;% if more channels than this, don't show electrode locations by default
EMARKER = '.';          % mark electrode locations with small disks
ECOLOR = [0 0 0];       % default electrode color = black
EMARKERSIZE = [];       % default depends on number of electrodes, set in code
EMARKERLINEWIDTH = 1;   % default edge linewidth for emarkers
EMARKERSIZE1CHAN = 20;  % default selected channel location marker size
EMARKERCOLOR1CHAN = 'red'; % selected channel location marker color
EMARKER2CHANS = [];      % mark subset of electrode locations with small disks
EMARKER2 = 'o';          % mark subset of electrode locations with small disks
EMARKER2COLOR = 'r';     % mark subset of electrode locations with small disks
EMARKERSIZE2 = 10;      % default selected channel location marker size
EMARKER2LINEWIDTH = 1;
EFSIZE = get(0,'DefaultAxesFontSize'); % use current default fontsize for electrode labels
HLINEWIDTH = 2;         % default linewidth for head, nose, ears
BLANKINGRINGWIDTH = .035;% width of the blanking ring
HEADRINGWIDTH    = .007;% width of the cartoon head ring
SHADING = 'flat';       % default 'shading': flat|interp
shrinkfactor = [];      % shrink mode (dprecated)
intrad       = [];      % default interpolation square is to outermost electrode (<=1.0)
plotrad      = [];      % plotting radius ([] = auto, based on outermost channel location)
headrad      = [];      % default plotting radius for cartoon head is 0.5
squeezefac = 1.0;
MINPLOTRAD = 0.15;      % can't make a topoplot with smaller plotrad (contours fail)
VERBOSE = 'off';
MASKSURF = 'off';
CONVHULL = 'off';       % dont mask outside the electrodes convex hull
DRAWAXIS = 'off';
PLOTDISK = 'off';
ContourVals = Values;
PMASKFLAG   = 0;
COLORARRAY  = { [1 0 0] [0.5 0 0] [0 0 0] };
COLORARRAY2 = { [gb 0] [gb 1/4] [gb 2/4] [gb 3/4] [gb 1] };
gb = [0 0];


[tmpeloc, labels, Th, Rd, indices] = readlocs(chanlocs);
Th = pi/180*Th;                              % convert degrees to radians
allchansind = 1:length(Th);
plotchans = indices;
[x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans = abs(plotchans);   % reverse indicated channel polarities
allchansind = allchansind(plotchans);
Th        = Th(plotchans);
Rd        = Rd(plotchans);
x         = x(plotchans);
y         = y(plotchans);
labels    = labels(plotchans); % remove labels for electrodes without locations
labels    = strvcat(labels); % make a label string matrix
if ~isempty(Values) && length(Values) > 1
    Values      = Values(plotchans);
    ContourVals = ContourVals(plotchans);
end

% Read plotting radius from chanlocs
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
default_intrad = 1;     % indicator for (no) specified intrad
intrad = min(1.0,max(Rd)*1.02);             % default: just outside the outermost electrode location

% set radius head cartoon
headrad = rmax;  % (anatomically correct)

% find plotting channels
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
intchans = find(Rd <= intrad); % interpolate channels in the radius intrad circle only

% % eliminate channels not plotted
allx      = x;
ally      = y;
intchans; % interpolate using only the 'intchans' channels
pltchans; % plot using only indicated 'plotchans' channels
allchansind = allchansind(pltchans);
intTh = Th(intchans);           % eliminate channels outside the interpolation area
intRd = Rd(intchans);
intx  = x(intchans);
inty  = y(intchans);
Th    = Th(pltchans);              % eliminate channels outside the plotting area
Rd    = Rd(pltchans);
x     = x(pltchans);
y     = y(pltchans);
labels= labels(pltchans,:);

% Squeeze channel locations to <= rmax
squeezefac = rmax/plotrad;
intRd = intRd*squeezefac; % squeeze electrode arc_lengths towards the vertex
Rd = Rd*squeezefac;       % squeeze electrode arc_lengths towards the vertex to plot all inside the head cartoon
intx = intx*squeezefac;   
inty = inty*squeezefac;  
x    = x*squeezefac;    
y    = y*squeezefac;   
allx    = allx*squeezefac;    
ally    = ally*squeezefac;   



% plot head outline
if headrad > 0                         % if cartoon head to be plotted
    headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
    heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];

    if ~ischar(HEADCOLOR) || ~strcmpi(HEADCOLOR,'none')
        %ringh= patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR,'linewidth', HLINEWIDTH); hold on
        headx = [rx(:)' rx(1)]*hin;
        heady = [ry(:)' ry(1)]*hin;
        ringh= plot(headx,heady);
        set(ringh, 'color',HEADCOLOR,'linewidth', HLINEWIDTH); hold on
    end

    % plot ears and nose
    base  = rmax-.0046;
    basex = 0.18*rmax;                   % nose width
    tip   = 1.15*rmax;
    tiphw = .04*rmax;                    % nose tip half width
    tipr  = .01*rmax;                    % nose tip rounding
    q = .04; % ear lengthening
    EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
    EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
    sf    = headrad/plotrad;                                          % squeeze the model ears and nose
    % by this factor
    if ~ischar(HEADCOLOR) || ~strcmpi(HEADCOLOR,'none')
        plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
            2*ones(size([basex;tiphw;0;-tiphw;-basex])),...
            'Color',HEADCOLOR,'LineWidth',HLINEWIDTH);                 % plot nose
        plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)    % plot left ear
        plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)   % plot right ear
    end

end

% show electrode info
plotax = gca;
axis square                                           % make plotax square
axis off
pos = get(gca,'position');
xlm = get(gca,'xlim');
ylm = get(gca,'ylim');
% textax = axes('position',pos,'xlim',xlm,'ylim',ylm);  % make new axes so clicking numbers <-> labels
% will work inside head cartoon patch
% axes(textax);
axis square                                           % make textax square
pos = get(gca,'position');
set(plotax,'position',pos);
xlm = get(gca,'xlim');
set(plotax,'xlim',xlm);
ylm = get(gca,'ylim');
set(plotax,'ylim',ylm);                               % copy position and axis limits again
axis equal;
lim = [-0.525 0.525];
%lim = [-0.56 0.56];
set(gca, 'xlim', lim); set(plotax, 'xlim', lim);
set(gca, 'ylim', lim); set(plotax, 'ylim', lim);
set(gca, 'xlim', lim); set(plotax, 'xlim', lim);
set(gca, 'ylim', lim); set(plotax, 'ylim', lim);

% z value for plotting electrode information (above the surf)
ELECTRODE_HEIGHT = 2.1;  

% electrode locations only
if strcmp(ELECTRODES,'on')   % plot electrodes as spots
    if isempty(EMARKER2CHANS)
        hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
            EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    else % plot markers for normal chans and EMARKER2CHANS separately
        hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
            EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
        hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
            EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
    end

    % electrode labels only
elseif strcmp(ELECTRODES,'labels')  % print electrode names (labels)
    for i = 1:size(labels,1)
        text(double(y(i)),double(x(i)),...
            ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','center',...
        	'VerticalAlignment','middle','Color',ECOLOR,...
        	'FontSize',EFSIZE)
    end

    % electrode locations + labels
elseif strcmp(ELECTRODES,'labelpoint')
    if isempty(EMARKER2CHANS)
        hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
            EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    else
        hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
            EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
        hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
            EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
    end
    for i = 1:size(labels,1)
        hh(i) = text(double(y(i)+0.01),double(x(i)),...
            ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','left',...
        	'VerticalAlignment','middle','Color', ECOLOR,'userdata', num2str(allchansind(i)), ...
        	'FontSize',EFSIZE, 'buttondownfcn', ...
    	    ['tmpstr = get(gco, ''userdata'');'...
   	     'set(gco, ''userdata'', get(gco, ''string''));' ...
         'set(gco, ''string'', tmpstr); clear tmpstr;'] );
    end
    
    % electrode locations + numbers
elseif strcmp(ELECTRODES,'numpoint')
    if isempty(EMARKER2CHANS)
        hp2 = plot3(y,x,ones(size(x))*ELECTRODE_HEIGHT,...
            EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    else
        hp2 = plot3(y(mark1chans),x(mark1chans),ones(size((mark1chans)))*ELECTRODE_HEIGHT,...
            EMARKER,'Color',ECOLOR,'markersize',EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
        hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
            EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
    end
    for i = 1:size(labels,1)
        hh(i) = text(double(y(i)+0.01),double(x(i)),...
            ELECTRODE_HEIGHT,num2str(allchansind(i)),'HorizontalAlignment','left',...
        	'VerticalAlignment','middle','Color', ECOLOR,'userdata', labels(i,:) , ...
        	'FontSize',EFSIZE, 'buttondownfcn', ...
    	    ['tmpstr = get(gco, ''userdata'');'...
   	     'set(gco, ''userdata'', get(gco, ''string''));' ...
         'set(gco, ''string'', tmpstr); clear tmpstr;'] );
    end
   
    % electrode numbers only
elseif strcmp(ELECTRODES,'numbers')
    for i = 1:size(labels,1)
        text(double(y(i)),double(x(i)),ELECTRODE_HEIGHT,int2str(allchansind(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','Color',ECOLOR,'FontSize',EFSIZE)
    end

    % emarker2 electrodes only  
elseif strcmp(ELECTRODES,'off') && ~isempty(EMARKER2CHANS)
    hp2b = plot3(y(mark2chans),x(mark2chans),ones(size((mark2chans)))*ELECTRODE_HEIGHT,...
        EMARKER2,'Color',EMARKER2COLOR,'markerfacecolor',EMARKER2COLOR,'linewidth',EMARKER2LINEWIDTH,'markersize',EMARKERSIZE2);
end

% Mark specified electrode locations with red filled disks
if strcmpi(STYLE,'blank') % if mark-selected-channel-locations mode
    for kk = 1:length(1:length(x))
        if abs(Values(kk))
            if strcmpi(PLOTDISK, 'off')
                angleRatio = real(Values(kk))/(real(Values(kk))+imag(Values(kk)))*360;
                radius     = real(Values(kk))+imag(Values(kk));
                allradius  = [0.02 0.03 0.037 0.044 0.05];
                radius     = allradius(radius);
                hp2 = disk(y(kk),x(kk),radius, [1 0 0], 0 , angleRatio, 16);
                if angleRatio ~= 360
                    hp2 = disk(y(kk),x(kk),radius, [0 0 1], angleRatio, 360, 16);
                end
            else
                tmpcolor = COLORARRAY{max(1,min(Values(kk), length(COLORARRAY)))};
                hp2 = plot3(y(kk),x(kk),ELECTRODE_HEIGHT,EMARKER,'Color', tmpcolor, 'markersize', EMARKERSIZE1CHAN);
                hp2 = disk(y(kk),x(kk),real(Values(kk))+imag(Values(kk)), tmpcolor, 0, 360, 10);
            end
        end
    end
end
