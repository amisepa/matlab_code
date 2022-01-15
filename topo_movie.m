    %% MY TOPO MOVIE
    % figure; pop_spectopo(EEG, 1, [], 'EEG', 'freq', [6 10 22], 'freqrange', [1 60],'electrodes','on');
    % pop_eegplot(EEG, 1, 1, 1);
    
    wind = EEG.srate;     %window size
    overlap = 50;           %(in %)
    Fs = EEG.srate;         %sample rate (in Hz)
    freqRange = 1:30;           %frequency of interest to compute
    
    % EEG = eeg_regepochs(EEG);
    
    % [psd, freqs] = get_psd(EEG.data,EEG.srate, overlap, Fs, freqRange, 'psd');
    
    disp('computing spectra...');
    for iWindow = 1:wind:size(EEG.data,2)
        [psd(:,:,iWindow), freqs] = get_psd(EEG.data(:,:,iWindow), EEG.srate, overlap, Fs, freqRange, 'psd');
    end
    
    figure;
    for iWindow = 1.5:0.5:10
        figure(2); colorbar;
        
        %     s1 = subplot(2,1,1);
        %     [psd(:,:,iEpoch), freqs] = get_psd(EEG.data(:,:,iEpoch), EEG.srate, overlap, Fs, freqRange, 'psd');
        %     for iChan = 1:EEG.nbchan
        %         psd_smooth(iChan,:,iEpoch) = conv(psd(iChan,:,iEpoch), ones(1,5)/5, 'same');
        %     end
        %     topodata = psd(:,:,iEpoch) - mean(psd(:,:,iEpoch),2);
        %     alpha = mean(topodata(:,8:13),2);
        delta = mean(psd(:,1:2,iWindow-1:iWindow),2);
        topoplot(delta(1:EEG.nbchan),EEG.chanlocs,'maplimits','absmax','numcontour', 0);
        hold on;
        title(['Delta power at ' num2str(iEpoch) ' s']);
        
        %     %asymmetry lines
        %     s2 = subplot(2,1,2);
        %     yline(0);
        % %     asy_front(:,iEpoch) = alpha(3)-alpha(2);
        % %     asy_temp(:,iEpoch) = alpha(4)-alpha(1);
        %     asy_front(:,iEpoch) = delta(3)-delta(2);
        %     asy_temp(:,iEpoch) = delta(4)-delta(1);
        %     xAxis = 1:size(EEG.data,3);
        %     plot(xAxis(iEpoch-1:iEpoch),asy_front(:,iEpoch-1:iEpoch),'b','LineWidth',1)
        %     hold on;
        %     plot(xAxis(iEpoch-1:iEpoch),asy_temp(:,iEpoch-1:iEpoch),'r','LineWidth',1)
        %     title(['Asymmetry levels at ' num2str(iEpoch) ' s']); legend('frontal','temporal');
        %
        %     text(double(xAxis(iEpoch)), double(asy_front(:,iEpoch)), sprintf('%s', num2str(round(asy_front(:,iEpoch)))), 'Parent', s2);
        %     text(double(xAxis(iEpoch)), double(asy_temp(:,iEpoch)), sprintf('%s', num2str(round(asy_temp(:,iEpoch)))), 'Parent', s2);
        
        %     pause(0.1)
    end
    hold off;
    
    mean(asy_front)
    mean(asy_temp)

    %% Simple 2-D movie
    % eeglab; close; % add path
    % eeglabp = fileparts(which('eeglab.m'));
    % EEG = pop_loadset(fullfile(eeglabp, 'sample_data', 'eeglab_data_epochs_ica.set'));
    
    % EEG = eeg_regepochs(EEG);
    
    % Above, convert latencies in ms to data point indices
    % pnts1 = round(eeg_lat2point(1, 1, EEG.srate, [EEG.xmin EEG.xmax]));
    % pnts2 = round(eeg_lat2point(256, 1, EEG.srate, [EEG.xmin EEG.xmax]));
    % scalpERP = mean(EEG.data(:,pnts1:pnts2),3);
    scalpERP = mean(EEG.data,3);
    
    % Smooth data
    for iChan = 1:size(scalpERP,1)
        scalpERP(iChan,:) = conv(scalpERP(iChan,:), ones(1,5)/5, 'same');
    end
    
    % 2-D movie
    figure;
    % [Movie,Colormap] = eegmovie(scalpERP, EEG.srate, EEG.chanlocs, 'framenum', 'off', 'vert', 0, 'startsec', -0.1, 'topoplotopt', {'numcontour' 0});
    [Movie,Colormap] = eegmovie(EEG.data, EEG.srate, EEG.chanlocs, 'framenum', 'off', 'vert', 0, 'startsec', -0.1, 'topoplotopt', {'numcontour' 0});
    seemovie(Movie,-5,Colormap);
    
    % save movie
    % vidObj = VideoWriter('erpmovie2d.mp4', 'MPEG-4');
    % open(vidObj);
    % writeVideo(vidObj, Movie);
    % close(vidObj);
    
    %% Simple 3-D movie
    % Use the graphic interface to coregister your head model with your electrode positions
    headplotparams1 = { 'meshfile', 'mheadnew.mat'       , 'transform', [0.664455     -3.39403     -14.2521  -0.00241453     0.015519     -1.55584           11      10.1455           12] };
    headplotparams2 = { 'meshfile', 'colin27headmesh.mat', 'transform', [0          -13            0          0.1            0        -1.57         11.7         12.5           12] };
    headplotparams  = headplotparams1; % switch here between 1 and 2
    
    % set up the spline file
    headplot('setup', EEG.chanlocs, 'STUDY_headplot.spl', headplotparams{:}); close
    
    % check scalp topo and head topo
    figure; headplot(scalpERP(:,end-50), 'STUDY_headplot.spl', headplotparams{:}, 'maplimits', 'absmax', 'lighting', 'on');
    figure; topoplot(scalpERP(:,end-50), EEG.chanlocs);
    figure('color', 'w'); [Movie,Colormap] = eegmovie( scalpERP, EEG.srate, EEG.chanlocs, 'framenum', 'off', 'vert', 0, 'startsec', -0.1, 'mode', '3d', 'headplotopt', { headplotparams{:}, 'material', 'metal'}, 'camerapath', [-127 2 30 0]);
    seemovie(Movie,-5,Colormap);
    
    % save movie
    vidObj = VideoWriter('erpmovie3d1.mp4', 'MPEG-4');
    open(vidObj);
    writeVideo(vidObj, Movie);
    close(vidObj);
    
    %% Using topoplot to make movie frames
    vidObj = VideoWriter('erpmovietopoplot.mp4', 'MPEG-4');
    open(vidObj);
    counter = 0;
    for latency = -100:10:600 %-100 ms to 1000 ms with 10 time steps
        figure; pop_topoplot(EEG,1,latency, 'My movie', [] ,'electrodes', 'off'); % plot'
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
        close;  % close current figure
    end
    close(vidObj);
