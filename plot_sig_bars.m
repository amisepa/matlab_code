%input h: binary vector or matrix of the same dimensionality as pvalues.  
%If the ith element of h is 1, then the ith p-value of pvalues is significant.  
%If the ith element of h is 0, then the ith p-value of pvalues is NOT significant.
%Cedric Cannard, May 2021

function plot_sig_bars(h)

count = 1;
for i = 1:length(h)
    if h(i) == 1 && i == 1
        seg(count,1) = 1;
    elseif h(i-1) == 0 && h(i) == 1
        seg(count,1) = i;
    end
    if h(i) == 1 && i == length(h)
        seg(count,2) = i; break        
    end
    if h(i) == 1 && h(i+1) == 0
        seg(count,2) = i;
        count = count+1;
    end
end

for i = 1:length(seg)
    hold on; plot([seg(i,1) seg(i,2)],[-20 -20],'k','LineWidth',3);
end

end