function plotmusedata(ax, data, srate, xl, yl)

colors = { 'r' [0 0.5 0] 'b' 'k' };
axes(ax); 
p1 = plot(linspace(xl(1), xl(2), srate*(xl(2)-xl(1))), data(1, xl(1)*srate+1:srate*xl(2))'); hold on;
p2 = plot(linspace(xl(1), xl(2), srate*(xl(2)-xl(1))), data(2, xl(1)*srate+1:srate*xl(2))'); 
p3 = plot(linspace(xl(1), xl(2), srate*(xl(2)-xl(1))), data(3, xl(1)*srate+1:srate*xl(2))'); 
p4 = plot(linspace(xl(1), xl(2), srate*(xl(2)-xl(1))), data(4, xl(1)*srate+1:srate*xl(2))');

set(p1, 'color', colors{1});
set(p2, 'color', colors{2});
set(p3, 'color', colors{3});
set(p4, 'color', colors{4});

ylim(yl); 
set(gca, 'yticklabel', []); 
xlim(xl);
