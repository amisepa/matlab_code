%% Converts raw EEG signals to Fourier domain to remove DC offset by 
% simple de-meaning signal (better than frequency domain). 
% 
% Usage:  EEG = rm_dc_offset(EEG);
% 
% Cedric Cannard, 2021

function EEG = rm_dc_offset(EEG)

disp('Removing DC drifts...');
% for iChan = 1:EEG.nbchan
%     ft = fft(EEG.data(iChan,:));
%     ft(1) = 0;                      % zero out the DC component
%     EEG.data(iChan,:) = ifft(ft);   % Inverse-transform back to time domain.
%     if mean(real(EEG.data(iChan,:))) > .005 || mean(real(EEG.data(iChan,:))) < -.005
%         errordlg('Mean should be closer to 0, DC drift removal must have failed.');
%         return
%     end
% end

% EEG.data = bsxfun(@minus, EEG.data, trimmean(EEG.data,20,2));
EEG.data = bsxfun(@minus, EEG.data, mean(EEG.data,2));
