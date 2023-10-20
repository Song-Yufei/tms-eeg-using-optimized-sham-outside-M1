% This script illustrates the steps of cleaning pipeline in the paper

% example subject 031 has 9 (3x3) datasets: three sessions x three targets (sma, mpfc, ag)
% For each target, the raw continuous EEG was first epoched around stimulus (-1.5s 1.5s), 
% and the sham events (150 trials) and active events(150 trials) were concatenated(Merged) in eeglab before the following cleaning steps.
% required toolbox, functions and files:
% fieldtrip(fieldtrip-20230613);
% eeglab(v2023.0): TESA(v1.2.0),FastICA_25 should be installed in the plugins;
% requiredFunction folder;
% original_chanlocs.mat
%% Step 1: Import data to EEGLAB
% load data
 eeglab;
 EEG = pop_loadset('filepath','\example_dataset_before_cleaning\',...
'filename', 'CHODME_031_3_SMA_TEPs(PRE)_Merge_epochs.set');

[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
 EEG = eeg_checkset(EEG);
 eeglab redraw     
%% Step 2: Remove unused channels 
EEG = pop_select( EEG, 'nochannel',{'APBr', 'FDIr','VEOG','HEOG'});
%% Step 3: Set baseline correction 
EEG = pop_rmbase( EEG, [-1000 -5]);
%% Step 4: Detrend (can take some time) and re-epoch (-1 to 1s) to aovid edge artifact following detrending
polynomial_order = 3; 
warning('off') % disable warnings 
EEG =  robust_detrend_EEG(EEG, polynomial_order);
warning('on') 
EEG = pop_epoch( EEG, {  }, [-1 1], 'epochinfo', 'yes');
%% Step 5: Decay Removal using Exponential fit (can take some time)
% define some parameters
timeToZP = [-4 14]./1000; % index2 corresponds to the starting point of decay search
baselineCorrTime = [-800 -10]./1000; % in second
time = EEG.times./1000; 
timeZero = find(time>=0,1);
tolTime = 3e-3; 
timeToFit = [];
noChangeTimeStart = 0.5; % from 500 msec. it should not change anymore 
weightTime = 40e-3; % time window or samples following tcut2 that model should prioritize, in s 
DecayFitSignal = EEG.data;

% Loop through trials and channels
 ft_progress('init', 'text', 'UnDecaying in progress...');
for trialNum=1:EEG.trials
     ft_progress(trialNum/EEG.trials, '\nUnDecaying is applied on all channels in trial %d from %d\n', trialNum, EEG.trials);            
    for chanNum=1:EEG.nbchan
       chanInd = chanNum;
       trialInd = trialNum;
       signal = double(EEG.data(chanInd,:,trialInd));
       [cleanSignalUD, ~, ~, ~] = unDecaying(signal, timeToZP, baselineCorrTime, ...
           EEG.srate, time, timeZero, timeToFit, weightTime, tolTime,noChangeTimeStart, [],[]);
       DecayFitSignal(chanNum,:,trialNum) = cleanSignalUD;
    end
end

EEG.data = DecayFitSignal;
EEG.timeToZP = timeToZP*1000;
EEG.baselineCorrTime = baselineCorrTime*1000;
%% Step 6: Remove TMS pulse artifact and interpolate before downsample
EEG.artCut = [EEG.timeToZP(1) EEG.timeToZP(2)+ EEG.tolTime];
EEG = pop_tesa_removedata( EEG, EEG.artCut );
EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );
%% Step 7: Downsample data (5000 Hz to 1000 Hz)
EEG = pop_resample( EEG, 1000);
%% Step 8: Remove heavily contaminated trial and channels (Visual inspection)
% ####Bad Channels#### 
% visualize
pop_eegplot( EEG, 1, 1, 1);
% ---insert channel labels to remove
EEG.Badchan = {'P5','P6','AF7'}; % these three channels are removed from all subjects' data due a technical error
EEG = pop_select( EEG, 'nochannel',EEG.Badchan);

% ####Bad Trials####
% visualize
EEG = pop_jointprob(EEG,1,1:size(EEG.data,1) ,5,3,0,0);
pop_rejmenu(EEG,1);
pause_script = input('Highlight bad trials, update marks and then press enter');
% in the GUI, select scroll data in the top right corner. 
% hit update marks and close the main trial rejection window. 
% hit enter in the command window.

EEG.BadTr = unique([find(EEG.reject.rejjp==1) find(EEG.reject.rejmanual==1)]);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 
%% Step 9: ICA without demean (only identify save ocular artifact topographies) 
EEG = pop_tesa_removedata( EEG, EEG.artCut);
EEG = runICA_TEP_on_EEG(EEG, []); % use defaut 35 ICs
% plot topographies and write down the indx for eye blink and movement.
figure; 
for i=1:35
    subplot(5,7,i)
    topoplot(EEG.icawinv(:,i), EEG.chanlocs);
end

EEG = pop_editset(EEG, 'setname', ['_ica.set']);
[ALLEEG,EEG,~] = eeg_store(ALLEEG,EEG); 
EEG = eeg_checkset(EEG); 
eeglab redraw 
%% Step 10: SOUND-SSP joint approach
% This SOUND-SSP joint approach estimates the signal subspace containing the TMS-related artifacts 
% and suppresses them from EEG signals. A framework can be referred in the
% paper: Hernandez-Pavon, J. C., Kugiumtzis, D., Zrenner, C., Kimiskidis, V. K., & Metsomaa, J. (2022). 
% Removing artifacts from TMS-evoked EEG: A methods review and a unifying theoretical framework.
%% Step 10.1: Define some parameters
clearvars data0 Pexp2 Uexp2 Utot data_corr corrM data_sound_filt data_sound Msound A

data0 = EEG.data;
times = EEG.times;
baselineCorrTime = EEG.baselineCorrTime;
eegchanlocs = EEG.chanlocs;
% artifacts time for defining artifacts subspace in ms (dataset dependent)
TimeArtifacts = [-10,90]; %in ms
TimeArtifacts_In = dsearchn(times',TimeArtifacts(1));
TimeArtifacts_End = dsearchn(times',TimeArtifacts(2));
% time windows to plot
TimePlots = [0,250]; %in ms
TimePlots_In = dsearchn(times',TimePlots(1));
TimePlots_End = dsearchn(times',TimePlots(2));
% computing lead-field matrix using spherical head model
LFM_sphere = ComputeSphericalLFM_chanlocs(EEG.chanlocs, 'CPz'); % CPz is the original recording reference
%% Step 10.2: Artifact subspace estimation 
% laplacian to maximize artifact-signal ratio, svd to decompose 
data0 = data0 - mean(data0(:,dsearchn(times',baselineCorrTime(1)):dsearchn(times',baselineCorrTime(2)),:),2); 
Xin=mean(data0(:,:,1:end),3);
Yrs=bnp_detrendLaplace(double(Xin(:,:,:)),5);

[Pexp2, Uexp2]=ProjectionMatrixComputation(Yrs(:,TimeArtifacts_In:TimeArtifacts_End,:), 1:EEG.pnts,[],...
    [TimePlots_In TimePlots_End], eegchanlocs, Xin);
% 1st row of plot: number of components to be chosen (Variance); raw data before any PCs to be deleted
% 2nd row of plot: topoplot of different PCs
% 3rd row of plot: the average time courses of different PCs
% 4th line of plot: all channel traces after 1:k PCs have been removed from the data. 

% left click mouse to toggle the visualization
% right click mouse to choose the first top k-dimensions with svd

%% 
% ---Inset 1:k to represent artifact subspace. large k can cause real signal suppressed
Utot=double([Uexp2(:,1:4)]); 
[Utot,~,~]=svds(Utot,size(Utot,2));  
%% Step 10.3: SOUND-SSP
% performs simultaneous SSP_SIR and SOUND: SSP is performed withing SOUND iterations without mixing the channels
data=data0; 
lambda =.005; 
iter = 10;
tsnoise = dsearchn(times',EEG.artCut(2)): TimeArtifacts_End; 
[data_corr, sigmas,dn, corrM] = SOUND_fast_SSP(data(:,tsnoise,:), LFM_sphere, iter,lambda, Utot,data(:,:,:));
% Re-reference to the averaged channels
data_corr = data_corr - mean(data_corr,1); 
data_sound = data_corr; 
% Re-reference the correction matrix
Msound = (eye(size(data_sound,1))-ones(size(data_sound,1),1)*ones(size(data_sound,1),1)'/size(data_sound,1))*corrM; 
%% Step 11: Remove ocular artifacts using beamforming
% ---insert ICA component indx for ocular artifacts from step 9
EEG.ocularICA = [2 9];
A = EEG.icawinv(:, EEG.ocularICA); 
% select 90% variance
[data_sound_filt, ~]=FilterOutGivenPatternsSound_BF(A, data_sound, Msound);
%% Step 12: Interpolate removed channels and convert to fieldtrip structure
% ---load original_chanlocs.mat
EEG = pop_editset(EEG, 'setname', ['_cleaned.set'] );
[ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG);
EEG = eeg_checkset(EEG);
eeglab redraw   

EEG.data = [];
EEG.data = data_sound_filt;
% interpolate missing data around TMS pulse artifacts
EEG = pop_tesa_removedata( EEG, [-5 EEG.artCut(2)] ); 
EEG = pop_tesa_interpdata( EEG, 'cubic', [5,5] );
% interpolate misssing chanels
EEG = pop_interp( EEG,original_chanlocs, 'spherical'); 
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 
EEG = eeg_checkset(EEG);
eeglab redraw     

% change event name, required by eeglab2ft function
EEG = rename_events(EEG, {'A - Stimulation'}, '1'); % active
EEG = rename_events(EEG, {'B - Stimulation'}, '2'); % sham
data_SIR_FT = eeglab2ft(EEG); 
save([filepath filename],'data_SIR_FT','-v7.3') 
                        