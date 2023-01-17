%% Pipieline to process continuous EEG data and generate a study
%
%  Cedric Cannard, 2022

clear; close all; clc;
dataDir = 'G:\Shared drives\Grants\Post Award Grants\(736) Bial Full-trance 2017\Research\Data\EEG\trance_bids_2021';
load(fullfile('G:\Shared drives\Grants\Post Award Grants\(736) Bial Full-trance 2017\Research\Interconnectivity\code\cedric_sInfo.mat'))
outDir = 'Desktop\trance';
mkdir(outDir)

cd(dataDir)
tmpSub = dir;
tmpSub = {tmpSub(9:end-1).name}';
cd(outDir)
eeglab; close;
% pop_editoptions('option_parallel', 0);
pop_editoptions('option_single', 0); % ensure double precision

counter1 = 1;
commands = {};
for iSub = 1%:5%13
    for iSes = 1:2

        disp('--------------------------------------------')
        disp(['Processing subject ' num2str(iSub) ' session ' num2str(iSes)])
        disp('--------------------------------------------')

        % file names and paths to load
        filepath = fullfile(dataDir, sprintf('sub-%3.3d',iSub), sprintf('ses-%2.2d',iSes), 'eeg');
        cd(filepath)
        tmp = dir; tmp = {tmp.name}';
        filename = tmp(contains(tmp, '.set'));
        
        % file names and path to save
        newpath = fullfile(outDir, sprintf('sub-%2.2d',iSub));
        if ~exist('newpath','dir'), mkdir(newpath); end
        newname = sprintf('sub-%2.2d.set',iSub);

        sInfo2(counter1).id = iSub;
        sInfo2(counter1).session = iSes;
        
        % import, channel locations, highpass filter
        EEG = pop_loadset('filename',filename,'filepath',filepath);
%         EEG = eeg_eegrej( EEG, [20 91679;320310 2247680]);     % for tutorial

        % Channel locations
%         chanlocpath = fileparts(which('dipfitdefs.m'));
%         chanlocfile = fullfile(chanlocpath,'standard_BEM','elec','standard_1005.elc');
        chanlocpath = fileparts(which('csd_transform.m'));
        chanlocfile = fullfile(chanlocpath, 'chanlocs_standard_BEM_1005.ced');
        EEG = pop_chanedit(EEG,'rplurchanloc',1,'lookup',chanlocfile);

        % Downsample and highpass filter
        EEG = pop_resample(EEG,256);
        order = firwsord('hamming',EEG.srate,0.5);     % get filter order for transition bandwidth = 0.5 Hz
        EEG = pop_firws(EEG,'forder',order,'fcutoff',.75,'ftype','highpass', ...
             'wtype','hamming','minphase',false,'plotfresp',false);
        
        % Average reference for ASR and ICA
        EEG = pop_reref(EEG,[]);
        
        % CSD-transformation
%         EEG = csd_transform(EEG,chanlocfile);

        % Detect flat channels
        maxFlat = 10;   % max flat segment tolerated in s (default = 5)
        maxJitter = 20; % max jitter tolerated during flatlines (as a multiple of epsilon; default - 20)
        oriEEG = EEG;
        badchannels = false(1,EEG.nbchan);
        for ichan = 1:EEG.nbchan
            zero_intervals = reshape(find(diff([false abs(diff(EEG.data(ichan,:)))<(maxJitter*eps) false])),2,[])';
            if max(zero_intervals(:,2) - zero_intervals(:,1)) > maxFlat*EEG.srate
                badchannels(ichan) = true; 
            end
        end
        
        %  Using clean_channels (correlation to its robust random estimate)
        minCorr = 0.8;      % minimum correlation between channels (default = .85)
        lineThresh = 5;     % line noise threshold (default = 4)
        winLength = 10;     % length of windows (in s) to compute corrThresh (default = 5)
        brokenTime = 0.33;  % max time (fraction of recording) of broken channel (0.1-0.6)
        nSamples = 100;     % ransac samples to generate random sampling consensus (in s; default = 50; higher is more robust but longer)
        [cleanEEG, rmchans] = clean_channels(EEG,minCorr,lineThresh,winLength,brokenTime,nSamples);
%     badChans = { oriEEG.chanlocs(rmchans).labels };
%     if EEG.nbchan < 59
%         error('More than 5 channels removed!')
%     end
%     vis_artifacts(EEG,oriEEG);
% %     idx = ~contains({oriEEG.chanlocs.labels}, {EEG.chanlocs.labels});
        sInfo2(counter1).badChans = {oriEEG.chanlocs(idx).labels};
%         toc(t)
        EEG = eeg_checkset(EEG);
%         EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical'); % interpolate
        figure; topoplot([],EEG.chanlocs,'style','blank', 'electrodes','labelpoint','chaninfo',EEG.chaninfo);
        saveas(gcf,fullfile(newpath, [sprintf('sub-%2.2d',iSub) '_ses-' num2str(iSes) '_bad_channels.png'])); close(gcf)

        % Remove large artifacts
        cutoff = 40;          % variance threshold (2-80)
        reconstruct = false;
        useriemannian = false;
        usegpu = false;
        m = memory;
        maxmem = round(.9 * (m.MemAvailableAllArrays/1000000),1);  % use half of available memory in MB
        disp(['Using 90% of available memory (' num2str(round(maxmem/1000,1)) ' GB)'])
        cleanEEG = clean_asr(EEG,cutoff,[],[],[],[],[],[],usegpu,useriemannian,maxmem);
%         [EEG,~,cleanEEG] = clean_artifacts(EEG,'Highpass','off','ChannelCriterion','off','LineNoiseCriterion','off', ...
%               'BurstRejection',~reconstruct,'BurstCriterion',cutoff, 'MaxMem',maxmem);

        % Remove segments from data 
        if reconstruct
            EEG = cleanEEG;
        else
            % segments of data that have changed
            mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
            EEG.etc.clean_sample_mask = ~mask;
           
            % get latency bounds of segments
            badData = reshape(find(diff([false mask false])),2,[])';
            badData(:,2) = badData(:,2)-1;
            
            % check if there are relevant events inside these bad portions
            tags = [];
            evLats = [EEG.event.latency];
            for iSeg = 1:size(badData,1)
                tags = find(ismember(evLats, badData(iSeg,1):badData(iSeg,2)));
                 if ~isempty(tags)
                     for iTag = 1:length(tags)
                        EEG.event(tags(iTag)).latency = badData(iSeg,1)-1;  % change latency to 1 sample before bad portion to preserve it
                        warning(['event ' num2str(tags(iTag)) ' was found in bad portion of data and moved at onset of bad portion to be preserved!']);
                     end
                 end
            end
            EEG = eeg_checkset(EEG,'eventconsistency');
    
            % remove segments
            EEG = pop_select(EEG,'nopoint',badData);
            fprintf('%3.1f %% of data were removed. \n', EEG.xmax/oriEEG.xmax*100)

        end
        EEG = eeg_checkset(EEG);

        % Remove irrecoverable time windows based on power
        disp('Now doing final post-cleanup of the output.');
        EEG = clean_windows(EEG,.25,[-inf 7]); 
        EEG = eeg_checkset(EEG);
        vis_artifacts(EEG,oriEEG); title(filename);
        
        % Interpolate bad channels
        EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical'); % interpolate

        % ICA
        datarank = sum(eig(cov(double(EEG.data'))) > 1E-7);
%         EEG = pop_runica(EEG, 'icatype','runica', 'extended',1,'pca',datarank);
        EEG = pop_runica(EEG,'icatype','picard','pca',datarank);
        EEG = pop_iclabel(EEG,'default');
        EEG = pop_icflag(EEG,[NaN NaN;0.95 1;0.95 1;0.99 1;0.95 1;0.99 1;NaN NaN]);
        pop_selectcomps(EEG, 1:dataRank);
        saveas(gcf,fullfile(newpath, [sprintf('sub-%2.2d',iSub) '_ses-' num2str(iSes) '_ICA.png'])); %close(gcf)
        bad_ic = find(EEG.reject.gcompreject);      % tag bad components
        if ~isempty(bad_ic)
            EEG = pop_subcomp(EEG, bad_ic);          % remove bad components from data
        end
        
        % Low pass filter
        order = firwsord('hamming',EEG.srate,2);     % get filter order for transition bandwidth = 0.5 Hz
        EEG = pop_firws(EEG,'forder',order,'fcutoff',45,'ftype','lowpass', ...
             'wtype','hamming','minphase',false,'plotfresp',false);
        
        % remove boundary events
        if sum(contains({EEG.event.type}, 'boundary'),'all') > 0
            idx = strcmpi({EEG.event.type}, 'boundary');
            EEG = pop_editeventvals(EEG,'delete',find(idx));
        end

        % Remove end markers, transitions, and first 10 s of each trial
        counter2 = 1;
        for iEv = 1:length(EEG.event)
            % remove period before first trial
            if iEv == 1
                rmdata(counter2,:) = [1 EEG.event(iEv).latency-1];
                counter2 = counter2 + 1;
            end
            % remove first 10 s of each trial
            if contains(EEG.event(iEv).type, 'start')
                rmdata(counter2,:) = [ EEG.event(iEv).latency+1 EEG.event(iEv).latency+EEG.srate*10 ];
                counter2 = counter2 + 1;
            end
            % remove end markers and transition periods
            if contains(EEG.event(iEv).type, 'end')
                if iEv ~= length(EEG.event)
                    rmdata(counter2,:) = [ EEG.event(iEv).latency-1 EEG.event(iEv+1).latency-1 ];
                    counter2 = counter2 + 1;
                else
                    rmdata(counter2,:) = [ EEG.event(iEv).latency-1 EEG.pnts ];
                    counter2 = counter2 + 1;
                end
            end
        end
        EEG = eeg_eegrej(EEG, rmdata);
        EEG = eeg_checkset(EEG);
%         pop_eegplot(EEG,1,1,1)

        % remove boundaries, and check length of each trial
%         idx = contains({EEG.event.type}, 'boundary');
%         EEG.event(idx) = []; 
        trialength = [];
        EEG = eeg_checkset(EEG);
        for iEv = 2:length(EEG.event)
            trialength(iEv,:) = ( EEG.event(iEv).latency - EEG.event(iEv-1).latency ) / EEG.srate / 60;
            if trialength(iEv,:) < 3.5
                error(['Trial ' num2str(iEv) ' is ' num2str(trialength(iEv,:)) ' min long and should be at least 3.5 min!'])
            end
        end
        
        % rename markers to have the same for each condition
        idx = contains({EEG.event.type}, 'rest');
        for iEv = 1:length(EEG.event)
            if idx(iEv)
                EEG.event(iEv).type = 'rest';
            else
                EEG.event(iEv).type = 'trance';
            end
        end
        
        % Filter out line noise and 16 Hz artifact (from GSR)
%         figure('Color','white');
%         for iChan = 1:EEG.nbchan
%             disp(num2str(iChan))
%             [power, f] = get_psd(EEG.data(iChan,:),EEG.srate*2,'hamming',50,[],EEG.srate,[1 80],'psd');
%             plot(f,power); hold on;
%         end
%         EEG = pop_zapline_plus(EEG, 'noisefreqs','line','coarseFreqDetectPowerDiff',4, ...
%             'chunkLength',0,'adaptiveNremove',1,'fixedNremove',1);
%         EEG = pop_cleanline(EEG,'bandwidth',2,'chanlist',1:EEG.nbchan,'computepower',true, ...
%             'linefreqs',16,'newversion',false,'normSpectrum',false,'p',0.001,'pad',2,'plotfigures',false, ...
%             'scanforlines',false,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',true,'winsize',4,'winstep',1);

        % Visualize data and save figure
        figure('Color','white');
        f = []; power = []; %power2 = [];
        for iChan = 1:EEG.nbchan
            disp(num2str(iChan))
%             subplot(8,8,iChan)
            [power(iChan,:), f] = get_psd(EEG.data(iChan,:),EEG.srate*2,'hamming',50,[],EEG.srate,[1 100],'psd');
%             power2(iChan,:) = get_psd(EEG2.data(iChan,:),EEG2.srate*2,'hamming',50,[],EEG2.srate,[1 30],'psd');
            plot(f,power(iChan,:)); hold on; %plot(f,power2(iChan,:));
%             title(num2str(iChan)); %legend('average', 'csd'); 
            title(filename)

            alpha(iChan,:) = mean(power(iChan,f>=8 & f<=13),'omitnan');
            theta(iChan,:) = mean(power(iChan,f>=3 & f<=7),'omitnan');
            beta(iChan,:) = mean(power(iChan,f>=17 & f<=30),'omitnan');
%             alpha2(iChan,:) = mean(power2(iChan,f>=8 & f<=13),'omitnan');
        end
        saveas(gcf,fullfile(newpath, [sprintf('sub-%2.2d',iSub) '_ses-' num2str(iSes) '_PSD.png'])); close(gcf)
        figure('Color','white');
%         subplot(2,1,1)
        subplot(3,1,1)
        topoplot(theta, EEG.chanlocs); title(filename); 
        clim([min(theta) max(theta)]); cbar; colormap('parula'); xlabel('theta')
        subplot(3,1,2)
        topoplot(alpha, EEG.chanlocs); title(filename); 
        clim([min(alpha) max(alpha)]); cbar; colormap('parula'); xlabel('alpha')
        subplot(3,1,3)
        topoplot(beta, EEG.chanlocs); title(filename); 
        clim([min(beta) max(beta)]); cbar; colormap('parula'); xlabel('beta')
%         subplot(2,1,2)
%         topoplot(alpha2, EEG2.chanlocs, 'gridscale', 50); xlabel('csd')
%         clim([min(alpha2) max(alpha2)]); cbar; colormap('parula'); 
        saveas(gcf,fullfile(newpath, [sprintf('sub-%2.2d',iSub) '_ses-' num2str(iSes) '_Topo.png'])); close(gcf)

        % copy of data to merge sessions
        if iSes == 1
            EEG1 = EEG;
        else
            EEG2 = EEG;
        end

        % history file
        history(iFile).subject = iSub;
        history(iFile).session = iSes;
        history(iFile).bad_channels = badChans;
        history(iFile).bad_data = []; % ???
        history(iFile).bad_comps = bad_ic;
        history(iFile).trial_length = trialength;
        iFile = iFile + 1;
    
    end

    % Merge across the two sessions and remove boundary
    EEG = eeg_emptyset;
    EEG = pop_mergeset(EEG1, EEG2);
    EEG = eeg_checkset(EEG);
    idx = contains({EEG.event.type}, 'boundary');
    EEG.event(idx) = [];
    EEG = eeg_checkset(EEG);
    if length(EEG.event) ~= 12
        error('wrong number of events after merging. There should be 6 rest and 6 trance markers'); 
    end

    % Compute PSD and plot
    figure('Color','white');
    f = []; power = [];
    for iChan = 1:EEG.nbchan
        disp(num2str(iChan))
        [power(iChan,:), f] = get_psd(EEG.data(iChan,:),EEG.srate*2,'hamming',50,[],EEG.srate,[1 100],'psd');
        plot(f,power(iChan,:)); hold on; title(filename)

        delta(iChan,:) = mean(power(iChan,f>=0 & f<=2.5),'omitnan');
        theta(iChan,:) = mean(power(iChan,f>=3 & f<=7),'omitnan');
        alpha(iChan,:) = mean(power(iChan,f>=8 & f<=13),'omitnan');
        beta(iChan,:) = mean(power(iChan,f>=17 & f<=30),'omitnan');
        gamma(iChan,:) = mean(power(iChan,f>=31 & f<=43),'omitnan');
    end
    saveas(gcf,fullfile(newpath, [sprintf('sub-%2.2d',iSub) '_PSD.png'])); close(gcf)
    
    % Scalp topography
    figure('Color','white');
    subplot(2,3,1)  % Delta
    topoplot(delta, EEG.chanlocs, 'emarker', {'.','k',5,1} );
    title('Delta','FontSize',10); 
    clim([min(delta) max(delta)]); 
    colorbar; colormap('parula'); %ylabel(c,'Power spectral difference (log)','FontSize',12)
    subplot(2,3,2)  % Theta
    topoplot(theta, EEG.chanlocs, 'emarker', {'.','k',5,1} );
    title('Theta','FontSize',10); 
    clim([min(theta) max(theta)]); 
    colorbar; colormap('parula'); %ylabel(c,'Power spectral difference (log)','FontSize',12)
    subplot(2,3,3)  % Alpha
    topoplot(alpha, EEG.chanlocs, 'emarker', {'.','k',5,1} );
    title('Alpha','FontSize',10); 
    clim([min(alpha) max(alpha)]); 
    colorbar; colormap('parula'); %ylabel(c,'Power spectral difference (log)','FontSize',12)
    subplot(2,3,4)  % Beta
    topoplot(beta, EEG.chanlocs, 'emarker', {'.','k',5,1} );
    title('Beta','FontSize',10); 
    clim([min(beta) max(beta)]); 
    colorbar; colormap('parula'); %ylabel(c,'Power spectral difference (log)','FontSize',12)
    subplot(2,3,5)  % Gamma
    topoplot(gamma, EEG.chanlocs, 'emarker', {'.','k',5,1} );
    title('Gamma','FontSize',10); 
    clim([min(gamma) max(gamma)]); 
    colorbar; colormap('parula'); %ylabel(c,'Power spectral difference (log)','FontSize',12)
    saveas(gcf,fullfile(newpath, [sprintf('sub-%2.2d',iSub) '_Topo.png'])); close(gcf)

    % epoch from -0.1 to 210 s (i.e., 3.5 minutes)
    EEG = pop_epoch(EEG,{},[-0.001  210],'newname','Merged and epoched','epochinfo','yes');
    EEG = eeg_checkset(EEG);
%     pop_eegplot(EEG,1,1,1)

%         % keep only mindwandering and trance data (exclude first 10 s of each run)
%         evIdx = contains({EEG.event.type}, 'rest');
%         rest = nan(11,2); trance = nan(11,2);
%         for iEv = 1:2:length(evIdx)
%             if evIdx(iEv) == 1 && evIdx(iEv+1) == 1
%                 rest(iEv,:) = [ EEG.event(iEv).latency + EEG.srate*10
%                     EEG.event(iEv+1).latency - EEG.srate*1 ];
%                 % rest(iEv,:) = [ EEG.event(iEv).latency-1
%                 % EEG.event(iEv+1).latency+1 ];
%             elseif evIdx(iEv) == 0 && evIdx(iEv+1) == 0
%                 trance(iEv,:) = [ EEG.event(iEv).latency + EEG.srate*10
%                     EEG.event(iEv+1).latency - EEG.srate*1 ];
%             else
%                 error('check events')
%             end
%         end
%         rest(isnan(rest(:,1)),:) = []; trance(isnan(trance(:,1)),:) = [];
%         if iSes == 1
%             restdata1 = pop_select(EEG, 'point', rest);       % rest data
%             trancedata1 = pop_select(EEG, 'point', trance);     % trance data
%         else
%             restdata2 = pop_select(EEG, 'point', rest);       % rest data
%             trancedata2 = pop_select(EEG, 'point', trance);     % trance data
%         end
%     end
% 
% balance channels across sessions (take file with least channels as ref)
%     if restdata1.nbchan ~= restdata2.nbchan
%         if restdata1.nbchan > restdata2.nbchan
%             idx = ~ismember({restdata1.chanlocs.labels}, {restdata2.chanlocs.labels});
%             restdata1 = pop_select(restdata1,'nochannel',{restdata1.chanlocs(idx).labels});
%             trancedata1 = pop_select(trancedata1,'nochannel',{trancedata1.chanlocs(idx).labels});
%         else
%             idx = ~ismember({restdata2.chanlocs.labels}, {restdata1.chanlocs.labels});
%             restdata2 = pop_select(restdata2,'nochannel',{restdata2.chanlocs(idx).labels});
%             trancedata2 = pop_select(trancedata2,'nochannel',{trancedata2.chanlocs(idx).labels});
%         end
%     end
%     restdata1 = eeg_checkset(restdata1);
%     restdata2 = eeg_checkset(restdata2);
%     trancedata1 = eeg_checkset(trancedata1);
%     trancedata2 = eeg_checkset(trancedata2);
% 
%         EEG.data = bsxfun(@rdivide, bsxfun(@minus, EEG.data, mean(EEG.data,2)), std(EEG.data, [],2));
%         EEG2.data = bsxfun(@rdivide, bsxfun(@minus, EEG2.data, mean(EEG2.data,2)), std(EEG2.data, [],2));
%         vis_artifacts(EEG,EEG2)
%         figure('Color','white'); 
%         count = 1;
%         for iChan = 1:4:64
%             subplot(4,4,count)
%             plot(EEG.data(iChan,1:256*5)); hold on; 
%             plot(EEG2.data(iChan,1:256*5)); 
%             legend('average', 'csd'); title(num2str(iChan));
%             count = count+1;
%         end
        
    % CREATE STUDY and SAVE
    EEG.subject = sprintf('sub-%2.2d',iSub);
%     EEG.condition = {'rest' 'trance'};
%     EEG.epochs = size(EEG.data,3);
    EEG.trials = size(EEG.data,3);
    EEG.run = [];
    EEG.saved = 'no';
    pop_saveset(EEG,'filepath',newpath,'filename',newname);
    commands = [ commands(:)' 'index' iSub 'load' fullfile(newpath, newname) ];
    [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','trance_eeg','commands', commands, ...
        'updatedat','on','savedat','off','rmclust','off');
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

    disp('--------------------------------------------')
    disp(['Subject ' num2str(iSub) ' done.'])
    disp('--------------------------------------------')

    counter2 = counter2+1;

end

save(fullfile(outDir, 'history.mat'), 'history')
[STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','trance_eeg.study','filepath',outDir);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

cd(outDir)
gong

%% study design

STUDY = std_makedesign(STUDY,ALLEEG,1,'name','STUDY.design 1','delfiles','off', 'defaultdesign','off', ...
    'variable1','type','values1',{'rest','trance'},'vartype1','categorical', ...
    'subjselect',{'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06', ...
    'sub-07','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13'});
[STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');

% STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','Design 1','delfiles','off','defaultdesign','off', ...
%     'variable1','condition','values1',{'rest','trance'},'vartype1','categorical', ...
%     'subjselect',{'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06', ...
%     'sub-07','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13'});
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');

% Precompute FFT/PSD
[STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on', ...
    'rmicacomps','off', 'interp','off','recompute','on','spec','on', ...
    'specparams',{'specmode','psd','logtrials','on','freqrange',[1 45]});
[STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');

%% Plot results with EEGLAB STUDY mode
clear; close all; clc;
outDir = 'C:\Users\IONSLAB\Desktop\channeling_matlab\data\data_processed2';
eeglab
[STUDY, ALLEEG] = pop_loadstudy('filename','trance_eeg.study','filepath',outDir);

STUDY = pop_statparams(STUDY,'condstats','on','method','perm','mcorrect','fdr','alpha',0.05);
STUDY = pop_specparams(STUDY, 'plotconditions','together','freqrange',[1 45] ,'averagechan','on');
[STUDY, specdata, specfreqs, pgroup, pcond] = std_specplot(STUDY,ALLEEG, ...
    'channels',{ALLEEG(1).urchanlocs.labels},'design',1,'noplot','on');
cond1 = cell2mat(specdata(1)); % condition 1
cond2 = cell2mat(specdata(2)); % condition 2
plotHDI(specfreqs',cond1,cond2,'Mean','dependent',.05,cell2mat(pcond)','rest','trance','Mean PSD + 95% HDI');
saveas(gcf,fullfile(outDir, 'results_all-freqs.png'));
% std_plotcurve(specfreqs',specdata,'plotconditions','together','plotstderr','on','figure','on');
% legend('rest', 'trance'); title('Mean PSD + 95% HDI');

% Topography for each band
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[0 2.5],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{ALLEEG(1).urchanlocs.labels},'design',1);
colormap('parula')
saveas(gcf,fullfile(outDir, 'topo_delta.png'));
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[3 7],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{ALLEEG(1).urchanlocs.labels},'design',1);
colormap('parula')
saveas(gcf,fullfile(outDir, 'topo_theta.png'));
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[8 13],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{ALLEEG(1).urchanlocs.labels},'design',1);
colormap('parula')
saveas(gcf,fullfile(outDir, 'topo_alpha.png'));
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[18 30],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{ALLEEG(1).urchanlocs.labels},'design',1);
colormap('parula')
saveas(gcf,fullfile(outDir, 'topo_beta.png'));
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[31 38],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{ALLEEG(1).urchanlocs.labels},'design',1);
colormap('parula')
saveas(gcf,fullfile(outDir, 'topo_gamma-low.png'));
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[63 78],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{ALLEEG(1).urchanlocs.labels},'design',1);
colormap('parula')
saveas(gcf,fullfile(outDir, 'topo_gamma-high.png'));

gong