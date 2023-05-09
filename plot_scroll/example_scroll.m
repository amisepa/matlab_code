
function superscroll_obj = example_scroll
    %example of superscroll.
    %returns superscroll object.
    
    %figure handle (click on figure to move scrollbar with arrows keys)
    fig = figure('position',[100 300 600 400]); 
    
    %make the reference axes (referenceAx) you want to scroll through.
    refax=axes('parent',fig);
    
   
    %put arbitrary plot on reference axes.
    X = unique(randi(300,[100,1]),'sorted');
    Y = randn(length(X),1);
    ploh=plot(X,Y,'-*r','linewidth',2,'parent',refax); clear X Y
    
    %arbitrarily set the reference axes limits to small portion of the
    %total X/Y data of the plot handle.
    set(refax,'xlim',[0 80],'ylim',[-1 2]);

    
    %make an X axis and Y axis scroll bar
    axeslist = {'X','Y'};
    superscroll_obj=superscroll(refax,axeslist);
    
    %build scrollbar, add keypress and resize functions
    autoscrollbar(superscroll_obj,ploh);
    
    %make one scroll bar yellow
    yellow(superscroll_obj.axdyn(1));
    
    
    title(refax,'pan and zoom using arrow keys, shift and control')
    
end