function plot_data2(data,varargin)


try
    p1 = varargin{1};
    p2 = varargin{2};
    p3 = varargin{3};
catch
end

DEFAULT_AXIS_COLOR = 'k';         % X-axis, Y-axis Color, text Color
SPACING_UNITS_STRING = '';        % '\muV' for microvolt optional units for g.spacingI Ex. uV
DEFAULT_GRID_SPACING = 1;         % Grid lines every n seconds

switch data
    case 'drawp' % Redraw EEG and change position

        % this test help to couple plot_data2 windows
        if exist('p3', 'var')
            figh = p3;
            figure(p3);
        else
            figh = gcf;                          % figure handle
        end

        if strcmp(get(figh,'tag'),'dialog')
            figh = get(figh,'UserData');
        end
%         ax0 = findobj('tag','backeeg','parent',figh); % axes handle
        ax1 = findobj('tag','eegaxis','parent',figh); % axes handle

        g = get(figh,'UserData');
        data = get(ax1,'UserData');
        ESpacing = findobj('tag','ESpacing','parent',figh);   % ui handle
        EPosition = findobj('tag','EPosition','parent',figh); % ui handle
        if ~isempty(EPosition) && ~isempty(ESpacing)
            if g.trialstag(1) == -1
                g.time    = str2double(get(EPosition,'string'));
            else
                g.time    = str2double(get(EPosition,'string'));
                g.time    = g.time - 1;
            end
            g.spacing = str2double(get(ESpacing,'string'));
        end

        if p1 == 1
            g.time = g.time-g.winlength*0.9;     % << subtract one window length
        elseif p1 == 2
            g.time = g.time-fastif(g.winlength>=5, round(g.winlength/5), g.winlength/5);             % < subtract one second
        elseif p1 == 3
            g.time = g.time+fastif(g.winlength>=5, round(g.winlength/5), g.winlength/5);             % > add one second
        elseif p1 == 4
            g.time = g.time+g.winlength*0.9;     % >> add one window length
        end

        if g.trialstag ~= -1 % time in second or in trials
            multiplier = g.trialstag;
        else
            multiplier = g.srate;
        end

        % Update edit box
        g.time = max(0,min(g.time,ceil((g.frames-1)/multiplier)-g.winlength));
        if g.trialstag(1) == -1
            set(EPosition,'string',num2str(g.time));
        else
            set(EPosition,'string',num2str(g.time+1));
        end
        set(figh, 'userdata', g);

        lowlim = round(g.time*multiplier+1);
        highlim = round(min((g.time+g.winlength)*multiplier+2,g.frames));

        % Plot data and update axes
        if ~isempty(g.data2)
            switch lower(g.submean) % subtract the mean ?
                case 'on'
                    meandata = mean(g.data2(:,lowlim:highlim)');
                    if any(isnan(meandata))
                        meandata = nan_mean(g.data2(:,lowlim:highlim)');
                    end
                otherwise, meandata = zeros(1,g.chans);
            end
        else
            switch lower(g.submean) % subtract the mean ?
                case 'on'
                    meandata = mean(data(:,lowlim:highlim)');
                    if any(isnan(meandata))
                        meandata = nan_mean(data(:,lowlim:highlim)');
                    end
                otherwise, meandata = zeros(1,g.chans);
            end
        end
        if strcmpi(g.plotdata2, 'on')
            hold on;
        else
            cla(ax1);
        end

        oldspacing = g.spacing;
        if g.envelope
            g.spacing = 0;
        end

        % plot channels whose "badchan" field is set to 1.
        % Bad channels are plotted first so that they appear behind the good
        % channels in the plot_data2 figure window.
        for i = 1:g.chans
            if strcmpi(g.plotdata2, 'on')
                tmpcolor = [ 1 0 0 ];
            else 
                tmpcolor = g.color{ mod(i-1,length(g.color))+1 } ;
            end

            if isfield(g, 'eloc_file') && isfield(g.eloc_file, 'badchan') && g.eloc_file(g.chans-i+1).badchan
                tmpcolor = [ .85 .85 .85 ];
                plot(ax1, data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
                    'color', tmpcolor, 'clipping','on')
                plot(ax1, 1,mean(data(g.chans-i+1,lowlim:highlim) -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing),2),'<r','MarkerFaceColor','r','MarkerSize',6);
            end

        end

        % plot good channels on top of bad channels (if g.eloc_file(i).badchan = 0... or there is no bad channel information)
        if strcmpi(g.plotdata2, 'off')
            tmpcolor = g.color{mod(g.chans-i,length(g.color))+1};
        end

        if (isfield(g, 'eloc_file') && isfield(g.eloc_file, 'badchan') && ~g.eloc_file(g.chans-i+1).badchan) || ...
                (~isfield(g, 'eloc_file')) || (~isfield(g.eloc_file, 'badchan'))
            plot(ax1, bsxfun(@plus, data(end:-1:1,lowlim:highlim), g.spacing*(1:g.chans)'-meandata(end:-1:1)')' + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
                'color', tmpcolor, 'clipping','on');
            % % This speeds things up below but there are
            % % drawback the plot handles should be passed as argument instead of using
            % % persistent var; also this does not allow plotting 2 sets of data
            %         persistent hplot;
            %         if isempty(hplot)
            %             hplot = plot(bsxfun(@plus, data(end:-1:1,lowlim:highlim), g.spacing*[1:g.chans]'-meandata(end:-1:1)')' + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), ...
            %                 'color', tmpcolor, 'clipping','on');
            %         else
            %             data2plot = bsxfun(@plus, data(end:-1:1,lowlim:highlim), g.spacing*[1:g.chans]'-meandata(end:-1:1)')' + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing);
            %             for iChan = 1:g.chans
            %                 set(hplot(iChan), 'ydata', data2plot(:,iChan)');
            %             end
            %         end
        end

        % draw selected channels
        if ~isempty(g.winrej) && size(g.winrej,2) > 2
            for tpmi = 1:size(g.winrej,1) % scan rows
                if (g.winrej(tpmi,1) >= lowlim && g.winrej(tpmi,1) <= highlim) || ...
                        (g.winrej(tpmi,2) >= lowlim && g.winrej(tpmi,2) <= highlim)
                    abscmin = max(1,round(g.winrej(tpmi,1)-lowlim));
                    abscmax = round(g.winrej(tpmi,2)-lowlim);
                    maxXlim = get(gca, 'xlim');
                    abscmax = min(abscmax, round(maxXlim(2)-1));
                    for i = 1:g.chans
                        if g.winrej(tpmi,g.chans-i+1+5)
                            plot(ax1, abscmin+1:abscmax+1,data(g.chans-i+1,abscmin+lowlim:abscmax+lowlim) ...
                                -meandata(g.chans-i+1)+i*g.spacing + (g.dispchans+1)*(oldspacing-g.spacing)/2 +g.elecoffset*(oldspacing-g.spacing), 'color','r','clipping','on')
                        end
                    end
                end
            end
        end
        g.spacing = oldspacing;
        set(ax1, 'Xlim',[1 g.winlength*multiplier],...
            'XTick',1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1);
        set(ax1, 'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'));

        % ordinates: even if all elec are plotted, some may be hidden
        set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );

        if g.children ~= 0
            if ~exist('p2', 'var')
                p2 =[];
            end
            plot_data2( 'drawp', p1, p2, g.children);
            figure(figh);
        end

        % draw second data if necessary
        if ~isempty(g.data2)
            tmpdata = data;
            set(ax1, 'userdata', g.data2);
            g.data2 = [];
            g.plotdata2 = 'on';
            set(figh, 'userdata', g);
            plot_data2('drawp', 0);
            g.plotdata2 = 'off';
            g.data2 = get(ax1, 'userdata');
            set(ax1, 'userdata', tmpdata);
            set(figh, 'userdata', g);
        else
            plot_data2('drawb'); % draw background first
        end

    case 'drawb' % Draw background 
        % Redraw EEG and change position

        ax0 = findobj('tag','backeeg','parent',gcf); % axes handle
        ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
        ylims=ylim(ax0);

        g = get(gcf,'UserData');  % Data (Note: this could also be global)

        % Plot data and update axes
        cla(ax0);
        hold(ax0, 'on');

        % plot rejected windows
        if g.trialstag ~= -1
            multiplier = g.trialstag;
        else
            multiplier = g.srate;
        end

        % draw rejection windows
        lowlim = round(g.time*multiplier+1);
        highlim = round(min((g.time+g.winlength)*multiplier+1));
%         displaymenu = findobj('tag','displaymenu','parent',gcf);
        if ~isempty(g.winrej) && g.winstatus
            if g.trialstag ~= -1 % epoched data
                indices = find((g.winrej(:,1)' >= lowlim & g.winrej(:,1)' <= highlim) | ...
                    (g.winrej(:,2)' >= lowlim & g.winrej(:,2)' <= highlim));
                if ~isempty(indices)
                    tmpwins1 = g.winrej(indices,1)';
                    tmpwins2 = g.winrej(indices,2)';
                    if size(g.winrej,2) > 2
                        tmpcols  = g.winrej(indices,3:5);
                    else 
                        tmpcols  = g.wincolor;
                    end
                    try
                        [cumul, indicescount] = histc(tmpwins1, (min(tmpwins1)-1):g.trialstag:max(tmpwins2));
                    catch 
                        [cumul, indicescount] = myhistc(tmpwins1, (min(tmpwins1)-1):g.trialstag:max(tmpwins2));
                    end
                    count = zeros(size(cumul));
                    %if ~isempty(find(cumul > 1)), find(cumul > 1), end
                    for tmpi = 1:length(tmpwins1)
                        poscumul = indicescount(tmpi);
                        heightbeg = count(poscumul)/cumul(poscumul);
                        heightend = heightbeg + 1/cumul(poscumul);
                        count(poscumul) = count(poscumul)+1;
                        h = patch(ax0, [tmpwins1(tmpi)-lowlim tmpwins2(tmpi)-lowlim ...
                            tmpwins2(tmpi)-lowlim tmpwins1(tmpi)-lowlim], ...
                            [heightbeg heightbeg heightend heightend], ...
                            tmpcols(tmpi,:));  % this argument is color
                        set(h, 'EdgeColor', get(h, 'facecolor'))
                    end
                end
            else
                event2plot1 = find ( g.winrej(:,1) >= lowlim & g.winrej(:,1) <= highlim );
                event2plot2 = find ( g.winrej(:,2) >= lowlim & g.winrej(:,2) <= highlim );
                event2plot3 = find ( g.winrej(:,1) <  lowlim & g.winrej(:,2) >  highlim );
                event2plot  = union_bc(union(event2plot1, event2plot2), event2plot3);

                for tpmi = event2plot(:)'
                    if size(g.winrej,2) > 2
                        tmpcols  = g.winrej(tpmi,3:5);
                    else 
                        tmpcols  = g.wincolor;
                    end
                    h = patch(ax0, [g.winrej(tpmi,1)-lowlim g.winrej(tpmi,2)-lowlim ...
                        g.winrej(tpmi,2)-lowlim g.winrej(tpmi,1)-lowlim], ...
                        [0 0 1 1], tmpcols);
                    set(h, 'EdgeColor', get(h, 'facecolor'))
                end
            end
        end

        % draw events if any
        if strcmpi(g.plotevent, 'on')

            MAXEVENTSTRING = g.maxeventstring;
            if MAXEVENTSTRING<0
                MAXEVENTSTRING = 0;
            elseif MAXEVENTSTRING>75
                MAXEVENTSTRING=75;
            end
            AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];

            % find event to plot
            event2plot    = find ( g.eventlatencies >=lowlim & g.eventlatencies <= highlim );
            if ~isempty(g.eventlatencyend)
                event2plot2 = find ( g.eventlatencyend >= lowlim & g.eventlatencyend <= highlim );
                event2plot3 = find ( g.eventlatencies  <  lowlim & g.eventlatencyend >  highlim );
                event2plot  = setdiff(union(event2plot, event2plot2), event2plot3);
            end
            for index = 1:length(event2plot)
                %Just repeat for the first one
                if index == 1
                    EVENTFONT = ' \fontsize{10} ';
                end

                % draw latency line
                tmplat = g.eventlatencies(event2plot(index))-lowlim-1;
%                 tmph = plot(ax0, [ tmplat tmplat ], ylims, 'color', g.eventcolors{ event2plot(index) }, ...
%                     'linestyle', g.eventstyle { event2plot(index) }, ...
%                     'linewidth', g.eventwidths( event2plot(index) ) );

                % schtefan: add Event types text above event latency line
                %             EVENTFONT = ' \fontsize{10} ';
                %             ylims=ylim;
                evntxt = strrep(num2str(g.events(event2plot(index)).type),'_','-');
                if length(evntxt)>MAXEVENTSTRING, evntxt = [ evntxt(1:MAXEVENTSTRING-1) '...' ]; end % truncate
                try
                    tmph2 = text(ax0, tmplat, ylims(2)-0.005, [EVENTFONT evntxt], ...
                        'color', g.eventcolors{ event2plot(index) }, ...
                        'horizontalalignment', 'left',...
                        'rotation',90);
                catch 
                end

                % draw duration is not 0
                if g.ploteventdur && ~isempty(g.eventlatencyend) ...
                        && g.eventwidths( event2plot(index) ) ~= 2.5 % do not plot length of boundary events
                    tmplatend = g.eventlatencyend(event2plot(index))-lowlim-1;
                    if tmplatend ~= 0
                        tmplim = ylims;
                        tmpcol = g.eventcolors{ event2plot(index) };
                        h = patch(ax0, [ tmplat tmplatend tmplatend tmplat ], ...
                            [ tmplim(1) tmplim(1) tmplim(2) tmplim(2) ], ...
                            tmpcol );  % this argument is color
                        set(h, 'EdgeColor', 'none')
                    end
                end
            end
        else % JavierLC
            MAXEVENTSTRING = 10; % default
            AXES_POSITION = [0.0964286 0.15 0.842 0.75-(MAXEVENTSTRING-5)/100];
        end

        if g.trialstag(1) ~= -1

            % plot trial limits
            tmptag = lowlim:highlim;
            tmpind = find(mod(tmptag-1, g.trialstag) == 0);
            for index = tmpind
                plot(ax0, [tmptag(index)-lowlim tmptag(index)-lowlim], [0 1], 'b--');
            end
            alltag = tmptag(tmpind);

            % compute Xticks
            tagnum = (alltag-1)/g.trialstag+1;
            set(ax0,'XTickLabel', tagnum,'YTickLabel', [],...
                'Xlim',[0 g.winlength*multiplier-1],...
                'XTick',alltag-lowlim+g.trialstag/2, 'YTick',[], 'tag','backeeg');

            tagpos  = [];
            tagtext = [];
            if ~isempty(alltag)
                alltag = [alltag(1)-g.trialstag alltag alltag(end)+g.trialstag]; % add border trial limits
            else
                alltag = [ floor(lowlim/g.trialstag)*g.trialstag ceil(highlim/g.trialstag)*g.trialstag ]+1;
            end

            nbdiv = 20/g.winlength; % approximative number of divisions
            divpossible = [ 100000./[1 2 4 5] 10000./[1 2 4 5] 1000./[1 2 4 5] 100./[1 2 4 5 10 20]]; % possible increments
            [tmp, indexdiv] = min(abs(nbdiv*divpossible-(g.limits(2)-g.limits(1)))); % closest possible increment
            incrementpoint = divpossible(indexdiv)/1000*g.srate;

            % tag zero below is an offset used to be sure that 0 is included
            % in the absicia of the data epochs
            if g.limits(2) < 0 
                tagzerooffset  = (g.limits(2)-g.limits(1))/1000*g.srate+1;
            else                
                tagzerooffset  = -g.limits(1)/1000*g.srate;
            end
            if tagzerooffset < 0, tagzerooffset = 0; end

            for i=1:length(alltag)-1
                if ~isempty(tagpos) && tagpos(end)-alltag(i)<2*incrementpoint/3
                    tagpos  = tagpos(1:end-1);
                end
                if ~isempty(g.freqlimits)
                    tagpos  = [ tagpos linspace(alltag(i),alltag(i+1)-1, nbdiv) ];
                else
                    if tagzerooffset ~= 0
                        tmptagpos = alltag(i)+tagzerooffset:-incrementpoint:alltag(i);
                    else
                        tmptagpos = [];
                    end
                    tagpos  = [ tagpos [tmptagpos(end:-1:2) alltag(i)+tagzerooffset:incrementpoint:(alltag(i+1)-1)]];
                end
            end

            % find corresponding epochs
            if ~g.isfreq
                tmplimit = g.limits;
                tpmorder = 1E-3;
            else
                tmplimit = g.freqlimits;
                tpmorder = 1;
            end
            tagtext = eeg_point2lat(tagpos, floor((tagpos)/g.trialstag)+1, g.srate, tmplimit,tpmorder);
            set(ax1,'XTickLabel', tagtext,'XTick', tagpos-lowlim+1 );
        else
            set(ax0,'XTickLabel', [],'YTickLabel', [],...
                'Xlim',[0 g.winlength*multiplier],...
                'XTick',[], 'YTick',[], 'tag','backeeg');

%             if g.isfreq
%                 set(ax1, 'XTickLabel', num2str((g.freqs(1):DEFAULT_GRID_SPACING:g.freqs(end))'),...
%                     'XTick',[1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1]);
%             else
                set(ax1,'XTickLabel', num2str((g.time:DEFAULT_GRID_SPACING:g.time+g.winlength)'),...
                    'XTick',1:multiplier*DEFAULT_GRID_SPACING:g.winlength*multiplier+1);
%             end

            set(ax1, 'Position', AXES_POSITION) % JavierLC
            set(ax0, 'Position', AXES_POSITION) % JavierLC
        end

        % ordinates: even if all elec are plotted, some may be hidden
        set(ax0, 'ylim',ylims );
        set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );

    case 'draws'
        % Redraw EEG and change scale

        ax1 = findobj('tag','eegaxis','parent',gcf);         % axes handle
        g = get(gcf,'UserData');
        data = get(ax1, 'userdata');
        ESpacing = findobj('tag','ESpacing','parent',gcf);   % ui handle
        EPosition = findobj('tag','EPosition','parent',gcf); % ui handle
        if g.trialstag(1) == -1
            g.time    = str2num(get(EPosition,'string'));
        else
            g.time    = str2num(get(EPosition,'string'))-1;
        end
        g.spacing = str2num(get(ESpacing,'string'));

        orgspacing= g.spacing;
        if p1 == 1
            g.spacing= g.spacing+ 0.1*orgspacing; % increase g.spacing(5%)
        elseif p1 == 2
            g.spacing= max(0,g.spacing-0.1*orgspacing); % decrease g.spacing(5%)
        end
        if round(g.spacing*100) == 0
            maxindex = min(10000, g.frames);
            g.spacing = 0.01*max(max(data(:,1:maxindex),[],2),[],1)-min(min(data(:,1:maxindex),[],2),[],1);  % Set g.spacingto max/min data
        end

        % update edit box
        set(ESpacing,'string',num2str(g.spacing,4))
        set(gcf, 'userdata', g);
        plot_data2('drawp', 0);
        set(ax1,'YLim',[0 (g.chans+1)*g.spacing],'YTick',[0:g.spacing:g.chans*g.spacing])
        set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );

        % update scaling eye (I) if it exists
        eyeaxes = findobj('tag','eyeaxes','parent',gcf);
        if ~isempty(eyeaxes)
            eyetext = findobj('type','text','parent',eyeaxes,'tag','thescalenum');
            set(eyetext,'string',num2str(g.spacing,4))
        end

        return;

    case 'window'  % change window size
        % get new window length with dialog box
        g = get(gcf,'UserData');
        result       = inputdlg2( { fastif(g.trialstag==-1,'New window length (s):', 'Number of epoch(s):') }, 'Change window length', 1,  { num2str(g.winlength) });
        if size(result,1) == 0 return; end

        g.winlength = eval(result{1});
        set(gcf, 'UserData', g);
        plot_data2('drawp',0);
        return;

    case 'winelec'  % change channel window size
        % get new window length with dialog box
        fig = gcf;
        g = get(gcf,'UserData');
        result = inputdlg2({ 'Number of channels to display:' } , 'Change number of channels to display', 1,  { num2str(g.dispchans) });
        if size(result,1) == 0 return; end

        g.dispchans = eval(result{1});
        if g.dispchans<0 || g.dispchans>g.chans
            g.dispchans =g.chans;
        end
        set(gcf, 'UserData', g);
        plot_data2('updateslider', fig);
        plot_data2('drawp',0);
        plot_data2('scaleeye', [], fig);
        return;

    case 'emaxstring'  % change events' string length  ;  JavierLC
        % get dialog box
        g = get(gcf,'UserData');
        result = inputdlg2({ 'Max events'' string length:' } , 'Change events'' string length to display', 1,  { num2str(g.maxeventstring) });
        if size(result,1) == 0 return; end
        g.maxeventstring = eval(result{1});
        set(gcf, 'UserData', g);
        plot_data2('drawb');
        return;

    case 'loadelect' % load channels
        [inputname,inputpath] = uigetfile('*','Channel locations file');
        if inputname == 0 return; end
        if ~exist([ inputpath inputname ])
            error('no such file');
        end

        AXH0 = findobj('tag','eegaxis','parent',gcf);
        plot_data2('setelect',[ inputpath inputname ],AXH0);
        return;

    case 'setelect'
        % Set channels
        eloc_file = p1;
        axeshand = p2;
        outvar1 = 1;
        if isempty(eloc_file)
            outvar1 = 0;
            return
        end

        tmplocs = readlocs(eloc_file);
        YLabels = { tmplocs.labels };
        YLabels = strvcat(YLabels);

        YLabels = flipud(char(YLabels,' '));
        set(axeshand,'YTickLabel',YLabels)

    case 'title'
        % Get new title
        h = findobj('tag', 'plot_data2title');

        if ~isempty(h)
            result       = inputdlg2( { 'New title:' }, 'Change title', 1,  { get(h(1), 'string') });
            if ~isempty(result), set(h, 'string', result{1}); end
        else
            result       = inputdlg2( { 'New title:' }, 'Change title', 1,  { '' });
            if ~isempty(result), h = textsc(result{1}, 'title'); set(h, 'tag', 'plot_data2title');end
        end

        return;

    case 'scaleeye'
        % Turn scale I on/off
        obj = p1;
        figh = p2;
        g = get(figh,'UserData');
        % figh = get(obj,'Parent');

        if ~isempty(obj)
            eyeaxes = findobj('tag','eyeaxes','parent',figh);
            children = get(eyeaxes,'children');
            if ischar(obj)
                if strcmp(obj, 'off')
                    set(children, 'visible', 'off');
                    set(eyeaxes, 'visible', 'off');
                    return;
                else
                    set(children, 'visible', 'on');
                    set(eyeaxes, 'visible', 'on');
                end
            else
                toggle = get(obj,'checked');
                if strcmp(toggle,'on')
                    set(children, 'visible', 'off');
                    set(eyeaxes, 'visible', 'off');
                    set(obj,'checked','off');
                    return;
                else
                    set(children, 'visible', 'on');
                    set(eyeaxes, 'visible', 'on');
                    set(obj,'checked','on');
                end
            end
        end

        eyeaxes = findobj('tag','eyeaxes','parent',figh);
        ax1 = findobj('tag','eegaxis','parent',gcf); % axes handle
        YLim = double(get(ax1, 'ylim'));

        ESpacing = findobj('tag','ESpacing','parent',figh);
        g.spacing= str2num(get(ESpacing,'string'));

        axes(eyeaxes); cla; axis off;
        set(eyeaxes, 'ylim', YLim);

        Xl = double([.35 .65; .5 .5; .35 .65]);
        Yl = double([ g.spacing g.spacing; g.spacing 0; 0 0] + YLim(1));
        plot(Xl(1,:),Yl(1,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline'); hold on;
        plot(Xl(2,:),Yl(2,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline');
        plot(Xl(3,:),Yl(3,:),'color',DEFAULT_AXIS_COLOR,'clipping','off', 'tag','eyeline');
        text(.5,(YLim(2)-YLim(1))/23+Yl(1),num2str(g.spacing,4),...
            'HorizontalAlignment','center','FontSize',10,...
            'tag','thescalenum')
        xlim([0 1]);
        if ~isempty(SPACING_UNITS_STRING)
            text(.5,-YLim(2)/23+Yl(4),SPACING_UNITS_STRING,...
                'HorizontalAlignment','center','FontSize',10, 'tag', 'thescale')
        end
        text(.5,(YLim(2)-YLim(1))/10+Yl(1),'Scale',...
            'HorizontalAlignment','center','FontSize',10, 'tag', 'thescale')
        set(eyeaxes, 'tag', 'eyeaxes');

    case 'noui'
        if ~isempty(varargin)
            plot_data2( varargin{:} ); fig = gcf;
        else
            fig = findobj('tag', 'PLOT_DATA2');
        end
        set(fig, 'menubar', 'figure');

        % find button and text
        obj = findobj(fig, 'style', 'pushbutton'); delete(obj);
        obj = findobj(fig, 'style', 'edit'); delete(obj);
        obj = findobj(fig, 'style', 'text');
        %objscale = findobj(obj, 'tag', 'thescale');
        %delete(setdiff(obj, objscale));
        obj = findobj(fig, 'tag', 'Eelec');delete(obj);
        obj = findobj(fig, 'tag', 'Etime');delete(obj);
        obj = findobj(fig, 'tag', 'Evalue');delete(obj);
        obj = findobj(fig, 'tag', 'Eelecname');delete(obj);
        obj = findobj(fig, 'tag', 'Etimename');delete(obj);
        obj = findobj(fig, 'tag', 'Evaluename');delete(obj);
        obj = findobj(fig, 'type', 'uimenu');delete(obj);
        set(gcf, 'paperpositionmode', 'manual', 'WindowButtonMotionFcn', []);

    case 'zoom' % if zoom
        fig = varargin{1};
        ax1 = findobj('tag','eegaxis','parent',fig);
        ax2 = findobj('tag','backeeg','parent',fig);
        tmpxlim  = get(ax1, 'xlim');
        tmpylim  = get(ax1, 'ylim');
        tmpxlim2 = get(ax2, 'xlim');
        set(ax2, 'xlim', get(ax1, 'xlim'));
        g = get(fig,'UserData');

        % deal with abscissa
        if g.trialstag ~= -1
            Eposition = str2num(get(findobj('tag','EPosition','parent',fig), 'string'));
            g.winlength = (tmpxlim(2) - tmpxlim(1))/g.trialstag;
            Eposition = Eposition + (tmpxlim(1) - tmpxlim2(1)-1)/g.trialstag;
            Eposition = round(Eposition*1000)/1000;
            set(findobj('tag','EPosition','parent',fig), 'string', num2str(Eposition));
        else
            Eposition = str2num(get(findobj('tag','EPosition','parent',fig), 'string'))-1;
            g.winlength = (tmpxlim(2) - tmpxlim(1))/g.srate;
            Eposition = Eposition + (tmpxlim(1) - tmpxlim2(1)-1)/g.srate;
            Eposition = round(Eposition*1000)/1000;
            set(findobj('tag','EPosition','parent',fig), 'string', num2str(Eposition+1));
        end

        % deal with ordinate
        g.elecoffset = tmpylim(1)/g.spacing;
        g.dispchans  = round(1000*(tmpylim(2)-tmpylim(1))/g.spacing)/1000;

        set(fig,'UserData', g);
        plot_data2('updateslider', fig);
        plot_data2('drawp', 0);
        plot_data2('scaleeye', [], fig);

        % reactivate zoom if 3 arguments
        if exist('p2', 'var') == 1
            if matVers < 8.4
                set(gcbf, 'windowbuttondownfcn', 'zoom(gcbf,''down''); plot_data2(''zoom'', gcbf, 1);');
            else
                % This is failing for us: http://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
                %               hManager = uigetmodemanager(gcbf);
                %               [hManager.WindowListenerHandles.Enabled] = deal(false);

                % Temporary fix
                wtemp = warning; warning off;
                set(gcbf, 'WindowButtonDownFcn', 'zoom(gcbf); plot_data2(''zoom'', gcbf, 1);');
                warning(wtemp);
            end
        end

    case 'updateslider' % if zoom
        fig = varargin{1};
        g = get(fig,'UserData');
        sliider = findobj('tag','eegslider','parent',fig);
        if g.elecoffset < 0
            g.elecoffset = 0;
        end
        if g.dispchans >= g.chans
            g.dispchans = g.chans;
            g.elecoffset = 0;
            set(sliider, 'visible', 'off');
        else
            set(sliider, 'visible', 'on');
            set(sliider, 'value', g.elecoffset/g.chans, ...
                'sliderstep', [1/(g.chans-g.dispchans) g.dispchans/(g.chans-g.dispchans)]);
            %'sliderstep', [1/(g.chans-1) g.dispchans/(g.chans-1)]);
        end
        if g.elecoffset < 0
            g.elecoffset = 0;
        end
        if g.elecoffset > g.chans-g.dispchans
            g.elecoffset = g.chans-g.dispchans;
        end
        set(fig,'UserData', g);
        plot_data2('scaleeye', [], fig);

    case 'drawlegend'
        fig = varargin{1};
        g = get(fig,'UserData');

        if ~isempty(g.events) % draw vertical colored lines for events, add event name text above
            nleg = length(g.eventtypes);
            fig2 = figure('numbertitle', 'off', 'name', '', 'visible', 'off', 'menubar', 'none', 'color', DEFAULT_FIG_COLOR);
            pos = get(fig2, 'position');
            set(fig2, 'position', [ pos(1) pos(2) 200 14*nleg+20]);

            for index = 1:nleg
                plot([10 30], [(index-0.5) * 10 (index-0.5) * 10], 'color', g.eventtypecolors{index}, 'linestyle', ...
                    g.eventtypestyle{ index }, 'linewidth', g.eventtypewidths( index )); hold on;
                if iscell(g.eventtypes)
                    th=text(35, (index-0.5)*10, g.eventtypes{index}, ...
                        'color', g.eventtypecolors{index}, 'interpreter', 'none');
                else
                    th=text(35, (index-0.5)*10, num2str(g.eventtypes(index)), ...
                        'color', g.eventtypecolors{index}, 'interpreter', 'none');
                end
            end
            xlim([0 130]);
            ylim([0 nleg*10]);
            axis off;
            set(fig2, 'visible', 'on');
        end


        % motion button: move windows or display current position (channel, g.time and activation)
        % case moved as subfunction
        % add topoplot
    case 'topoplot'
        fig = varargin{1};
        g = get(fig,'UserData');
        if ~isstruct(g.eloc_file) || ~isfield(g.eloc_file, 'theta') || isempty( [ g.eloc_file.theta ])
            return;
        end
        ax1 = findobj('tag','backeeg','parent',fig);
        tmppos = get(ax1, 'currentpoint');
        ax1 = findobj('tag','eegaxis','parent',fig); % axes handle
        
        % plot vertical line
        yl = ylim;
        plot([ tmppos tmppos ], yl, 'color', [0.8 0.8 0.8]);

        if g.trialstag ~= -1,
            lowlim = round(g.time*g.trialstag+1);
        else, lowlim = round(g.time*g.srate+1);
        end
        data = get(ax1,'UserData');
        datapos = max(1, round(tmppos(1)+lowlim));
        datapos = min(datapos, g.frames);

        figure; topoplot(data(:,datapos), g.eloc_file);
        if g.trialstag == -1
            latsec = (datapos-1)/g.srate;
            title(sprintf('Latency of %d seconds and %d milliseconds', floor(latsec), round(1000*(latsec-floor(latsec)))));
        else
            trial = ceil((datapos-1)/g.trialstag);

            latintrial = eeg_point2lat(datapos, trial, g.srate, g.limits, 0.001);
            title(sprintf('Latency of %d ms in trial %d', round(latintrial), trial));
        end
        return;

        % release button: check window consistency, add to trial boundaries
    case 'defupcom'
        fig = varargin{1};
        g = get(fig,'UserData');
        ax1 = findobj('tag','backeeg','parent',fig);
        g.incallback = 0;
        set(fig,'UserData', g);  % early save in case of bug in the following
        if strcmp(g.mocap,'on'), g.winrej = g.winrej(end,:);end % nima
        if ~isempty(g.winrej)', ...
                if g.winrej(end,1) == g.winrej(end,2) % remove unitary windows
                g.winrej = g.winrej(1:end-1,:);
                else
                    if g.winrej(end,1) > g.winrej(end,2) % reverse values if necessary
                        g.winrej(end, 1:2) = [g.winrej(end,2) g.winrej(end,1)];
                    end
                    g.winrej(end,1) = max(1, g.winrej(end,1));
                    g.winrej(end,2) = min(g.frames, g.winrej(end,2));
                    if g.trialstag == -1 % find nearest trials boundaries if necessary
                        I1 = find((g.winrej(end,1) >= g.winrej(1:end-1,1)) & (g.winrej(end,1) <= g.winrej(1:end-1,2)) );
                        if ~isempty(I1)
                            g.winrej(I1,2) = max(g.winrej(I1,2), g.winrej(end,2)); % extend epoch
                            g.winrej = g.winrej(1:end-1,:); % remove if empty match
                        else
                            I2 = find((g.winrej(end,2) >= g.winrej(1:end-1,1)) & (g.winrej(end,2) <= g.winrej(1:end-1,2)) );
                            if ~isempty(I2)
                                g.winrej(I2,1) = min(g.winrej(I2,1), g.winrej(end,1)); % extend epoch
                                g.winrej = g.winrej(1:end-1,:); % remove if empty match
                            else
                                I2 = find((g.winrej(end,1) <= g.winrej(1:end-1,1)) & (g.winrej(end,2) >= g.winrej(1:end-1,1)) );
                                if ~isempty(I2)
                                    g.winrej(I2,:) = []; % remove if empty match
                                end
                            end
                        end
                    end
                end
        end
        set(fig,'UserData', g);
        plot_data2('drawp', 0);
        if strcmp(g.mocap,'on'), show_mocap_for_plot_data2(g.winrej); g.winrej = g.winrej(end,:); end % nima

        % push button: create/remove window
    case 'defdowncom'
        fig = varargin{1};
        g = get(fig,'UserData');
        if strcmp(g.mocap,'on'), show_mocap_timer = timerfind('tag','mocapDisplayTimer'); if ~isempty(show_mocap_timer),  end; end % nima

        ax1 = findobj('tag','backeeg','parent',fig);
        tmppos = get(ax1, 'currentpoint');
        if strcmp(get(fig, 'SelectionType'),'normal')

            fig = varargin{1};
            g = get(fig,'UserData');
            ax1 = findobj('tag','backeeg','parent',fig);
            tmppos = get(ax1, 'currentpoint');
            g = get(fig,'UserData'); % get data of background image {g.trialstag g.winrej incallback}
            if g.incallback ~= 1 % interception of nestest calls
                if g.trialstag ~= -1
                    lowlim = round(g.time*g.trialstag+1);
                    highlim = round(g.winlength*g.trialstag);
                else
                    lowlim  = round(g.time*g.srate+1);
                    highlim = round(g.winlength*g.srate); % THIS IS NOT TRUE WHEN ZOOMING

                end
                if (tmppos(1) >= 0) && (tmppos(1) <= highlim)
                    if isempty(g.winrej) Allwin=0;
                    else Allwin = (g.winrej(:,1) < lowlim+tmppos(1)) & (g.winrej(:,2) > lowlim+tmppos(1));
                    end
                    if any(Allwin) % remove the mark or select electrode if necessary
                        lowlim = find(Allwin==1);
                        if g.setelectrode  % select electrode
                            ax2 = findobj('tag','eegaxis','parent',fig);
                            tmppos = get(ax2, 'currentpoint');
                            tmpelec = g.chans + 1 - round(tmppos(1,2) / g.spacing);
                            tmpelec = min(max(tmpelec, 1), g.chans);
                            g.winrej(lowlim,tmpelec+5) = ~g.winrej(lowlim,tmpelec+5); % set the electrode
                        else  % remove mark
                            g.winrej(lowlim,:) = [];
                        end
                    else
                        if g.trialstag ~= -1 % find nearest trials boundaries if epoched data
                            alltrialtag = 0:g.trialstag:g.frames;
                            I1 = find(alltrialtag < (tmppos(1)+lowlim) );
                            if ~isempty(I1) && I1(end) ~= length(alltrialtag)
                                g.winrej = [g.winrej' [alltrialtag(I1(end)) alltrialtag(I1(end)+1) g.wincolor zeros(1,g.chans)]']';
                            end
                        else
                            g.incallback = 1;  % set this variable for callback for continuous data
                            if size(g.winrej,2) < 5
                                g.winrej(:,3:5) = repmat(g.wincolor, [size(g.winrej,1) 1]);
                            end
                            if size(g.winrej,2) < 5+g.chans
                                g.winrej(:,6:(5+g.chans)) = zeros(size(g.winrej,1),g.chans);
                            end
                            g.winrej = [g.winrej' [tmppos(1)+lowlim tmppos(1)+lowlim g.wincolor zeros(1,g.chans)]']';
                        end
                    end
                    set(fig,'UserData', g);
                    % plot_data2('drawp', 0);  % redraw background
                end
            end
        end
    otherwise
        error(['Error - invalid plot_data2() parameter: ', data])
end

