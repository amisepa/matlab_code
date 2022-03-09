%Displays all channels and 2 min of data (30 s sections)
%Check bad channels to tag them as bad for rejection
% 
% Usage: badChan = tag_badChannels(EEG)
% 
% Cedric Cannard, 2020

function badChan = tag_badChannels(EEG)

disp(['File length: ' num2str(EEG.xmax)]);

if EEG.xmax < 61
    badChan = true(1,4);
    cprintf('blue', 'File too short, skipping to next file \n');
    return
end

badChan = false(1,4);

%Create axes
% if ~exist('s1', 'var')
% f = figure;
s1 = axes('unit', 'normalized', 'position' , [0.05 0.8 0.7 0.2]);
s2 = axes('unit', 'normalized', 'position' , [0.05 0.55 0.7 0.2]);
s3 = axes('unit', 'normalized', 'position' , [0.05 0.3 0.7 0.2]);
s4 = axes('unit', 'normalized', 'position' , [0.05 0.05 0.7 0.2]);
s5 = axes('unit', 'normalized', 'position' , [0.77 0.6 0.21 0.4]);
% end

% EEG = clean_drifts(EEG,[0.25 0.75]);   
% EEG = pop_eegfiltnew(EEG, 'hicutoff',45,'plotfreqz',0);
standDev = std(EEG.data, [], 2);
EEG.data = bsxfun(@rdivide, bsxfun(@minus, EEG.data, mean(EEG.data,2)), standDev*3);
EEG.data = bsxfun(@plus, EEG.data, [12 8 4 0]');

% Plot EEG.data in each axe:
plotmusedata(s1, EEG.data, EEG.srate, [0 30], [-2 14]); set(gca,'XTick',[]); xlabel('0 - 30 s'); %ylim([-standDev(2)/2 standDev(1)*2]);
plotmusedata(s2, EEG.data, EEG.srate, [30 60], [-2 14]); set(gca,'XTick',[]); xlabel('30 - 60 s'); %ylim([-standDev(2)/2 standDev(2)*2]);
plotmusedata(s3, EEG.data, EEG.srate, [60 90], [-2 14]); set(gca,'XTick',[]); xlabel('60 - 90 s'); %ylim([-standDev(3)/2 standDev(3)*2]);
if EEG.xmax >= 120
    plotmusedata(s4, EEG.data, EEG.srate, [90 120], [-2 14]); set(gca,'XTick',[]); xlabel('90 - 120 s'); %ylim([-standDev(4)/2 standDev(4)*2]);
end
axes(s5);spectopo(EEG.data, size(EEG.data,2), EEG.srate);xlim([.5 50]); ylabel(''); set(gca, 'yticklabel', []); %ylim([-60 20]);

% Create GUI text and buttons
uicontrol('unit', 'normalized', 'position', [0.77 0.41  0.2, 0.05], 'style', 'text', 'string', 'Select each bad channel you see (checked = bad channel):');
chan1 = uicontrol('unit', 'normalized', 'position', [0.78 0.35  0.2, 0.05], 'style', 'checkbox', 'string', sprintf('Channel 1  (black; SD: %1.0f)', standDev(1)));
chan2 = uicontrol('unit', 'normalized', 'position', [0.78 0.31 0.2, 0.05], 'style', 'checkbox', 'string', sprintf('Channel 2  (blue; SD: %1.0f)', standDev(2)));
chan3 = uicontrol('unit', 'normalized', 'position', [0.78 0.27 0.2, 0.05], 'style', 'checkbox', 'string', sprintf('Channel 3  (green; SD: %1.0f)', standDev(3)));
chan4 = uicontrol('unit', 'normalized', 'position', [0.78 0.23 0.2, 0.05], 'style', 'checkbox', 'string', sprintf('Channel 4  (red; SD: %1.0f)', standDev(4)));
nextButton = uicontrol('unit', 'normalized', 'position', [0.85 0.1 0.05, 0.05], 'style', 'pushbutton', 'string', 'Next', 'callback', 'set(gcbo, ''userdata'', 1);');
waitfor(nextButton, 'userdata');

% Get bad channels from check boxes and close figure
badChan = logical([chan1.Value chan2.Value chan3.Value chan4.Value]);
% badChan = {EEG.chanlocs(badChan).labels};
% close(f);

disp('Bad channels tagged.')

end