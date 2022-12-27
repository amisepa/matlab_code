%% Converts raw EEG signals to Fourier domain to remove DC component 
% (i.e., slow-frequency drifts) and converts back to time domain. 
% 
% Usage:    EEG = rm_DC(EEG);
% 
% Cedric Cannard, 2021

function EEG = rm_DC(EEG)

disp('Removing DC drifts...');
for iChan = 1:EEG.nbchan
    ft = fft(EEG.data(iChan,:));
    ft(1) = 0;                      % zero out the DC component
    EEG.data(iChan,:) = ifft(ft);   % Inverse-transform back to time domain.
    if mean(real(EEG.data(iChan,:))) > .005 || mean(real(EEG.data(iChan,:))) < -.005
        warning('Mean should be closer to 0, DC drift removal must have failed.')
    end
end
