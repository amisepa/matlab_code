%% COmpute functional connectivity

function [FC, freqs] = compute_fc(data, conn_meas)

% data = EEG.roi.source_roi_data;
% win_len = 2; % seconds
% conn_meas = {'CS'}; 
% fs = EEG.srate;
% nfft = 50; % number of frequency bins

nROI    = size(data,1);
% pnts    = size(data,2);
nTrials = size(data,3);

% number of PCs
nPCA = 3;
ninds = 0; inds = {};
for iroi = 1:nROI
    for jroi = (iroi+1):nROI
        inds{ninds+1} = {(iroi-1)*nPCA + 1:nPCA, (jroi-1)*nPCA + 1:nPCA};
        ninds = ninds + 1;
    end
end

% Connectivity - Coh/PDC/iCoh
for iTrial = 1:nTrials

    fprintf('Trial %g/%g\n',iTrial,nTrials)

    % Data for this epoch
    tmpdata = double(squeeze(data(:,:,iTrial)));

    % Weighted phase lag index (wPLI), Granger Causality (GC), Time-reversed Granger Causality (TRGC), 
    % Maximized imaginary coherency (MIC)
    if any(strcmpi(conn_meas,{'wPLI' 'GC' 'TRGC' 'MIC'}))
        fc = data2sctrgcmim(tmpdata,0,[],[],[],[],conn_meas);  
        if strcmpi(conn_meas, 'MIM') || strcmpi(conn_meas, 'MIC')
            MI = fc(:, :);
            fc = get_connect_mat( MI, nROI, +1);
        elseif  strcmpi(conn_meas, 'GC') || strcmpi(conn_meas, 'TRGC')
            TRGCnet = fc(:, :, 1) - fc(:, :, 2);
            fc = get_connect_mat( TRGCnet, nROI, -1); 
        else % wPLI
            warning(strcat("Only the first principal component will be used to determine ", conn_meas))
            measure = rm_components(fc, nPCA); % only keep the first principal component
            fc = measure;
        end
        FC(:,:,:,iTrial) = fc.conn_meas;
    
    % Cross spectrum (CS), Coherence (aCOH), Complex-valued coherence (cCOH), 
    % imaginary part of coherency (iCOH), Partial directed coherence (PDC), 
    % Directed transfer entropy (DTF), Time-reversed directed transfer
    % entropy (TRDTF), Phase-amplitude coupling (PAC)
    elseif any(strcmpi(conn_meas,{'CS' 'aCOH' 'cCOH' 'iCOH' 'PDC' 'DTF' 'TRDTF' 'PAC'}))
        [fc, freqs] = data2spwctrgc(tmpdata,[],[],[],[],inds,conn_meas);   
        % nlags = fc.nlags;
        fc = rm_components(fc, nPCA); % only keep the first principal component
        FC(:,:,:,iTrial) = abs(fc.(conn_meas{:}));
    else
        error("Unknown connectivity measure. Should be one from: 'CS' 'aCOH' 'cCOH' 'iCOH' 'PDC' 'DTF' 'TRDTF' 'PAC' 'wPLI' 'GC' 'TRGC' 'MIC'")
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%55 fdMVAR_5order %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Su = eye(nROI,nROI);
    % [dc,dtf,pdc,gpdc,coh,pcoh,pcoh2,h,s,pp,f] = fdMVAR_5order(tmpdata,Su,fs*2,fs);
    % coh = abs(coh);
    % idx = f >= 1 <= 3;

end
