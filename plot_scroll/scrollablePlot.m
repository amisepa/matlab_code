function [superscroll_obj,ploth] = scrollablePlot(plotProps,axeslist)
    %Use this function instead of plot() to make a plot with scrollbar(s)
    %e.g. scrollablePlot((1:4,rand(4,1),'parent',A,'color','k'),{'X','Y'})
    
    %axes list are the axes you want to make scroll bars for
    %{'X','Y'} makes X and Y scrollbar, {'X'} makes X scroll bar, {'Y'} makes Y scroll bar
    
    %plotProps are property value pairs for plot()
    
    %function returns superscroll_obj (superscroll object), and the handle
    %to plot you made with plotProps.
    
    ploth=plot(plotProps{:}); %plot handle
    ax=ploth.Parent;
    %set initial axes limits to a portion of the total plotted X Data
    %(unless the plot is empty or all NaN's)
    XX=ploth.XData; 
    if ~isempty(XX)
        YY=ploth.YData;
        if ~all(isnan(XX))
            totxlim=[nanmin(XX),nanmax(XX)];
            xlim=[totxlim(1),totxlim(1)+0.2*diff(totxlim)]; %arbitrarily show first 20%
        end
        %make all the Y data visible (unless Y data is all NaN's)
        if ~all(isnan(YY))
            ylim = [nanmin(YY),nanmax(YY)];
        end
    end
    set(ax,'xlim',xlim,'ylim',ylim);
    
    superscroll_obj = superscroll(ax,axeslist);
    %build scrollbars here 
    autoscrollbar(superscroll_obj,ploth)
end