%% Scan data to find bad trials using amplitude and high-frequency power 
%
% INPUT: 
%   EEG
%
% OUTPUT: 
%   badTrials
%
% Cedric Cannard, Dec 2022

function badTrials = find_badTrials(EEG)

fprintf('Looking for bad trials. This may take a few minutes... \n')
sigRMS = nan(size(EEG.data,3),1);
sigPower = nan(size(EEG.data,3),1);
progressbar('Looking for bad trials')
for iEpoch = 1:size(EEG.data,3)
    power = nan(EEG.nbchan,30);
    sigRMS(iEpoch,:) = rms(rms(EEG.data(:,:,iEpoch)));
    for iChan = 1:EEG.nbchan
        [power(iChan,:), f] = get_psd(EEG.data(iChan,:,iEpoch),EEG.srate,'hamming',50,[],EEG.srate,[70 100],'psd');
    end
    % hold on; plot(f,mean(power));
    sigPower(iEpoch,:) = rms(rms(power));  % mean across channels and RMS across frequencies

    progressbar(iEpoch/size(EEG.data,3))
end

% badAmp = find(sigRMS > 5*std(sigRMS));        % bad epochs based on amplitude
badAmp = find(isoutlier(sigRMS,'gesd'));        % 'mean' 'median' 'quartiles' 'grubbs' 'gesd' (default)
% badPow = find(sigPower < 5*std(sigPower));    % bad epochs based on high-frequency power
badPow = find(isoutlier(sigPower,'mean'));      % 'mean' (more lax) 'median' 'quartiles' 'grubbs' 'gesd'

badTrials = sort(unique([badAmp; badPow]));

gong
