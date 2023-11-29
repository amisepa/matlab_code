%% plot_timefreq
% Plot of time-frequency and ITC  
% 
% usage:
%   g = plot_timefreq(P, R, Pboot, Rboot, ERP, freqs, times, mbase, maskersp, maskitc, g);
% 
% Cedric Cannard, 2023
 
function plot_timefreq(pwr, itc, pwr_H0, itc_H0, freqs, times, mbase, maskersp, maskitc, params)

plotitc = true;
plotersp = true;
scale = 'log';
basenorm = 'off';
alpha = 0.05;
pcontour = 'off';

% auto scalling
% erspmax = [max(max(abs(P)))]/2;
% if strcmpi(scale,'abs') && strcmpi(basenorm,'off')
%     erspmax = [max(max(abs(pwr)))];
%     if erspmax > 1
%         erspmax = [1-(erspmax-1) erspmax];
%     else 
%         erspmax = [erspmax 1+(1-erspmax)];
%     end
% end

if plotitc

    if ~isreal(itc)
        Rangle = angle(itc);
        Rsign = sign(imag(itc));
        itc = abs(itc); % convert coherence vector to magnitude
        setylim = 1;
    else
        Rangle = zeros(size(itc)); % avoid error because Rangle does not exist
        Rsign = ones(size(itc));
        setylim = 0;
    end


    if plotersp
        ordinate1 = 0.67; 
        ordinate2 = 0.1; 
        height = 0.33;
    else
        ordinate2 = 0.1; 
        height = 0.9; 
    end
else
    ordinate1 = 0.1; 
    height = 0.9;
    if plotersp
        ordinate1 = 0.1; 
        height = 0.9;
    end
end

set(gcf,'DefaultAxesFontSize', 10)
pos = get(gca,'position');
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)];
axis off;

%% ERSP

if plotersp
    h(1) = axes('Position',[.1 ordinate1 .9 height].*s+q);
    set(h(1), 'tag', 'ersp');

    if strcmpi(scale, 'abs') && strcmpi(basenorm, 'off')
        baseval = 1;
    else 
        baseval = 0;
    end
    
    if ~isnan(alpha)
        if strcmpi(g.pcontour, 'off') && ~isempty(maskersp) % zero out nonsignif. power differences
            PP(~maskersp) = baseval;
            %PP = PP .* maskersp;
        elseif isempty(maskersp)
            if size(PP,1) == size(pwr_H0,1) && size(PP,2) == size(pwr_H0,2)
                PP(find(PP > pwr_H0(:,:,1) & (PP < pwr_H0(:,:,2)))) = baseval;
                pwr_H0 = squeeze(mean(pwr_H0,2));
                if size(pwr_H0,2) == 1, pwr_H0 = pwr_H0'; end
            else
                PP(find((PP > repmat(pwr_H0(:,1),[1 length(times)])) ...
                    & (PP < repmat(pwr_H0(:,2),[1 length(times)])))) = baseval;
            end
        end
    end
 
        % find color limits
        if isempty(g.erspmax)
            if g.ERSP_CAXIS_LIMIT == 0
                g.erspmax = [-1 1]*1.1*max(max(abs(pwr(:,:))));
            else
                g.erspmax = g.ERSP_CAXIS_LIMIT*[-1 1];
            end
        elseif length(g.erspmax) == 1
            g.erspmax = [ -g.erspmax g.erspmax];
        end
        if isnan( g.baseline(1) ) && g.erspmax(1) < 0
            g.erspmax = [ min(min(pwr(:,:))) max(max(pwr(:,:)))];
        end

        % plot image
        if ~strcmpi(g.freqscale, 'log')
            imagesc(times,freqs,PP(:,:), g.erspmax);
        else
            imagesclogy(times,freqs,PP(:,:),g.erspmax);
        end
        set(gca,'ydir',g.hzdir);  % make frequency ascend or descend

        % put contour for multiple comparison masking
        if ~isempty(maskersp) && strcmpi(g.pcontour, 'on')
            hold on; [tmpc tmph] = contour(times, freqs, maskersp);
            set(tmph, 'linecolor', 'k', 'linewidth', 0.25)
        end
        
        hold on
        plot([0 0],[0 freqs(end)],'--m','LineWidth',g.linewidth); % plot time 0
        if ~isnan(g.marktimes) % plot marked time
            for mt = g.marktimes(:)'
                plot([mt mt],[0 freqs(end)],'--k','LineWidth',g.linewidth);
            end
        end
        hold off
        set(h(1),'YTickLabel',[],'YTick',[])
        set(h(1),'XTickLabel',[],'XTick',[])
        if ~isempty(g.vert)
            for index = 1:length(g.vert)
                line([g.vert(index), g.vert(index)], [min(freqs) max(freqs)], 'linewidth', 1, 'color', 'm');
            end
        end

        h(2) = gca;
        h(3) = cbar('vert'); % ERSP colorbar axes
        set(h(2),'Position',[.1 ordinate1 .8 height].*s+q)
        set(h(3),'Position',[.95 ordinate1 .05 height].*s+q)
        title([ 'ERSP(' g.unitpower ')' ])

        %
        %%%%% plot marginal ERSP mean below ERSP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(4) = axes('Position',[.1 ordinate1-0.1 .8 .1].*s+q);

        E = [min(pwr(:,:),[],1);max(pwr(:,:),[],1)];

        % plotting limits
        if isempty(g.erspmarglim)
            g.erspmarglim = [min(E(1,:))-max(max(abs(E)))/3 max(E(2,:))+max(max(abs(E)))/3];
        end

        plot(times,E,[0 0],g.erspmarglim, '--m','LineWidth',g.linewidth)
        xlim([min(times) max(times)])
        ylim(g.erspmarglim)

        tick = get(h(4),'YTick');
        set(h(4),'YTick',[tick(1) ; tick(end)])
        set(h(4),'YAxisLocation','right')
        set(h(4),'TickLength',[0.020 0.025]);
        xlabel('Time (ms)')
        ylabel(g.unitpower)

        %
        %%%%% plot mean spectrum to left of ERSP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(5) = axes('Position',[0 ordinate1 .1 height].*s+q);

        if isnan(g.baseline)              % Ramon :for bug 1657
            E = zeros(size(freqs));
        else
            E = mbase;
        end

        if ~isnan(E(1))

            % plotting limits
            if isempty(g.speclim)
               % g.speclim = [min(E)-max(abs(E))/3 max(E)+max(abs(E))/3];
               if all(~isnan(mbase))
                   g.speclim = [min(mbase)-max(abs(mbase))/3 max(mbase)+max(abs(mbase))/3]; % RMC: Just for plotting
               else
                   g.speclim = [min(E)-max(abs(E))/3 max(E)+max(abs(E))/3];
               end
            end

            % plot curves
            if ~strcmpi(g.freqscale, 'log')
                plot(freqs,E,'LineWidth',g.linewidth); hold on;
                if ~isnan(g.alpha) && size(pwr_H0,2) == 2
                    try
                        plot(freqs,pwr_H0(:,:)'+[E;E], 'g', 'LineWidth',g.linewidth)
                        plot(freqs,pwr_H0(:,:)'+[E;E], 'k:','LineWidth',g.linewidth)
                    catch
                        plot(freqs,pwr_H0(:,:)+[E E], 'g', 'LineWidth',g.linewidth)
                        plot(freqs,pwr_H0(:,:)+[E E], 'k:','LineWidth',g.linewidth)
                    end
                end
                if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
                if g.speclim(1) ~= g.speclim(2), ylim(g.speclim); end; % Ramon :for bug 1657 

            else % 'log'
                semilogx(freqs,E,'LineWidth',g.linewidth); hold on;
                if ~isnan(g.alpha)
                    try
                        semilogx(freqs,pwr_H0(:,:)'+[E;E],'g', 'LineWidth',g.linewidth)
                        semilogx(freqs,pwr_H0(:,:)'+[E;E],'k:','LineWidth',g.linewidth)
                    catch
                        semilogx(freqs,pwr_H0(:,:)+[E E],'g', 'LineWidth',g.linewidth)
                        semilogx(freqs,pwr_H0(:,:)+[E E],'k:','LineWidth',g.linewidth)
                    end
                end
                if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
                if g.speclim(1) ~= g.speclim(2), ylim(g.speclim); end; %RMC
                set(h(5),'View',[90 90])
                divs = linspace(log(freqs(1)), log(freqs(end)), 10);
                set(gca, 'xtickmode', 'manual');
                divs = ceil(exp(divs)); divs = unique_bc(divs); % ceil is critical here, round might misalign
                set(gca, 'xtick', divs);
            end
            set(h(5),'TickLength',[0.020 0.025]);
            set(h(5),'View',[90 90])
            xlabel('Frequency (Hz)')
            if strcmp(g.hzdir,'normal')
                set(gca,'xdir','reverse');
            else
                set(gca,'xdir','normal');
            end
            ylabel(g.unitpower)
            tick = get(h(5),'YTick');
            if (length(tick)>2)
                set(h(5),'YTick',[tick(1) ; tick(end-1)])
            end
        end
end

switch lower(g.plotitc)
    case 'on'
        %
        %%%%%%%%%%%% Image the ITC %%%%%%%%%%%%%%%%%%
        %
        h(6) = axes('Position',[.1 ordinate2 .9 height].*s+q); % ITC image
        if ishandle(h(1));set(h(1), 'tag', 'itc');end

        if abs(itc(1,1)-1) < 0.0001, g.plotphaseonly = 'on'; end
        if strcmpi(g.plotphaseonly, 'on')
            RR = Rangle/pi*180;
        else
            RR = itc;
        end
        if ~isnan(g.alpha)
            if ~isempty(maskitc) && strcmpi(g.pcontour, 'off')
                RR = RR .* maskitc;
            elseif isempty(maskitc)
                if size(RR,1) == size(itc_H0,1) && size(RR,2) == size(itc_H0,2)
                    tmp = gcf;
                    if size(itc_H0,3) == 2	 RR(find(RR > itc_H0(:,:,1) & RR < itc_H0(:,:,2))) = 0;
                    else                   RR(find(RR < itc_H0)) = 0;
                    end
                    itc_H0 = mean(itc_H0(:,:,end),2);
                else
                    RR(find(RR < repmat(itc_H0(:),[1 length(times)]))) = 0;
                end
            end
        end

        if g.ITC_CAXIS_LIMIT == 0
            coh_caxis = min(max(max(itc(:,:))),1)*[-1 1]; % 1 WAS 0.4 !
        else
            coh_caxis = g.ITC_CAXIS_LIMIT*[-1 1];
        end

        if strcmpi(g.plotphaseonly, 'on')
            if ~strcmpi(g.freqscale, 'log')
                imagesc(times,freqs,RR(:,:)); % <---
            else
                imagesclogy(times,freqs,RR(:,:)); % <---
            end
            g.itcmax = [-180 180];
            setylim = 0;
        else
            if max(coh_caxis) == 0,              % toby 10.02.2006
                coh_caxis = [-1 1];
            end
            if ~strcmpi(g.freqscale, 'log')
                if exist('Rsign') && strcmp(g.plotphasesign, 'on')
                    imagesc(times,freqs,Rsign(:,:).*RR(:,:),coh_caxis); % <---
                else
                    imagesc(times,freqs,RR(:,:),coh_caxis); % <---
                end
            else
                if exist('Rsign') && strcmp(g.plotphasesign, 'on')
                    imagesclogy(times,freqs,Rsign(:,:).*RR(:,:),coh_caxis); % <---
                else
                    imagesclogy(times,freqs,RR(:,:),coh_caxis); % <---
                end
            end
        end
        set(gca,'ydir',g.hzdir);  % make frequency ascend or descend

        % plot contour if necessary
        if ~isempty(maskitc) && strcmpi(g.pcontour, 'on')
            hold on; [tmpc tmph] = contour(times, freqs, maskitc);
            set(tmph, 'linecolor', 'k', 'linewidth', 0.25)
        end

        if isempty(g.itcmax)
            g.itcmax = caxis;
        elseif length(g.itcmax) == 1
            g.itcmax = [ -g.itcmax g.itcmax ];
        end
        caxis(g.itcmax);

        hold on
        plot([0 0],[0 freqs(end)],'--m','LineWidth',g.linewidth);
        if ~isnan(g.marktimes)
            for mt = g.marktimes(:)'
                plot([mt mt],[0 freqs(end)],'--k','LineWidth',g.linewidth);
            end
        end
        hold off
        set(h(6),'YTickLabel',[],'YTick',[])
        set(h(6),'XTickLabel',[],'XTick',[])
        if ~isempty(g.vert)
            for index = 1:length(g.vert)
                line([g.vert(index), g.vert(index)], [min(freqs) max(freqs)], 'linewidth', 1, 'color', 'm');
            end
        end

        h(7) = gca;
        h(8) = cbar('vert');
        %h(9) = get(h(8),'Children'); % make the function crash
        set(h(7),'Position',[.1 ordinate2 .8 height].*s+q)
        set(h(8),'Position',[.95 ordinate2 .05 height].*s+q)
        if setylim
            set(h(8),'YLim',[0 g.itcmax(2)]);
        end
        if strcmpi(g.plotphaseonly, 'on')
            title('ITC phase')
        else
            title('ITC')
        end

        %
        %%%%% plot the ERP below the ITC image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %

        h(10) = axes('Position',[.1 ordinate2-0.1 .8 .1].*s+q); % ERP

        if isempty(g.erplim)
            ERPmax = max(ERP);
            ERPmin = min(ERP);
            g.erplim = [ ERPmin - 0.1*(ERPmax-ERPmin) ERPmax + 0.1*(ERPmax-ERPmin) ];
        end

        plot(ERPtimes,ERP, [0 0],g.erplim,'--m','LineWidth',g.linewidth);
        hold on;
        plot([times(1) times(length(times))],[0 0], 'k');
        xlim([min(ERPtimes) max(ERPtimes)]);
        ylim(g.erplim)
        set(gca,'ydir',g.ydir);

        tick = get(h(10),'YTick');
        set(h(10),'YTick',[tick(1) ; tick(end)])
        set(h(10),'TickLength',[0.02 0.025]);
        set(h(10),'YAxisLocation','right')
        xlabel('Time (ms)')
        ylabel('\muV')
        if (~isempty(g.topovec))
            if length(g.topovec) ~= 1, ylabel(''); end; % ICA component
        end
        E = nan_mean(itc(:,:)'); % don't let a few NaN's crash this

        %
        %%%%% plot the marginal mean left of the ITC image %%%%%%%%%%%%%%%%%%%%%
        %

        h(11) = axes('Position',[0 ordinate2 .1 height].*s+q); % plot the marginal mean
        % ITC left of the ITC image
        % set plotting limits
        if isempty(g.itcavglim)
            if ~isnan(g.alpha)
                g.itcavglim = [ min(E)-max(E)/3 max(itc_H0)+max(itc_H0)/3];
            else
                g.itcavglim = [ min(E)-max(E)/3 max(E)+max(E)/3];
            end
        end
        if max(g.itcavglim) == 0 || any(isnan(g.itcavglim))
            g.itcavglim = [-1 1];
        end
        
        % plot marginal ITC
        if ~strcmpi(g.freqscale, 'log')
            plot(freqs,E,'LineWidth',g.linewidth); hold on;
            if ~isnan(g.alpha)
                plot(freqs,itc_H0,'g', 'LineWidth',g.linewidth)
                plot(freqs,itc_H0,'k:','LineWidth',g.linewidth)
            end
            if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
            ylim(g.itcavglim)
        else
            semilogx(freqs,E,'LineWidth',g.linewidth); hold on;
            if ~isnan(g.alpha)
                semilogx(freqs,itc_H0(:),'g', 'LineWidth',g.linewidth)
                semilogx(freqs,itc_H0(:),'k:','LineWidth',g.linewidth)
            end
            if freqs(1) ~= freqs(end), xlim([freqs(1) freqs(end)]); end
            ylim(g.itcavglim)
            divs = linspace(log(freqs(1)), log(freqs(end)), 10);
            set(gca, 'xtickmode', 'manual');
            divs = ceil(exp(divs)); divs = unique_bc(divs); % ceil is critical here, round might misalign
            set(gca, 'xtick', divs);
         end

        % ITC plot details
        tick = get(h(11),'YTick');
        if length(tick) > 1
            set(h(11),'YTick',[tick(1) ; tick(length(tick))])
        end
        set(h(11),'View',[90 90])
        %set(h(11),'TickLength',[0.020 0.025]);
        xlabel('Frequency (Hz)')
        if strcmp(g.hzdir,'normal')
            set(gca,'xdir','reverse');
        else
            set(gca,'xdir','normal');
        end
        ylabel('ERP')

end %switch
