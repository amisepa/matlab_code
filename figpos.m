function figpos(h)

% PC's active screen size
screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

% Include figure toolbar and border
toolbar_height = 77;
window_border  = 5;

% The Format of Matlab is this:
%[left, bottom, width, height]
m_left   = pc_width/2 + 300 + window_border;
m_bottom = pc_height/2 + 100 - toolbar_height - 1;
m_height = .33*pc_height;
m_width  = .25*pc_width - 1;

%Set the correct position of the figure
set(h, 'Position', [m_left, m_bottom, m_width, m_height]);

% To print to file correctly, set these below and use the "-r72" scaling to 
% get the proper format.
% set(h, 'PaperUnits', 'points');
% set(h, 'PaperSize', [width, height]); 
% set(h, 'PaperPosition', [0, 0, width, height]); %[ left, bottom, width, height]
