% This script loads the preprocessed dataset after the tms-eeg cleaning pipeline, calculates stimulus-related averaging per subject, 
% and makes a grand average over all subjects 

% Each subject has 9 (3 x3) trial-level datasets, corresponding to three sessions (S1, S2, S3) 
% and three targets (SMA, AG, mPFC). Each dataset is merged with active ('1') 
% and sham ('2') events. The datasets are stored in the fieldtrip data structure. 
%%
clear, clc

subjies={'002','003', '004', '007', '009', '010', '011', '012',...
    '013', '014', '015', '017', '018', '019', '020', '021', '022','024','025','027','028','030','031','032'};

% ----load dataset of '' target
TARG='SMA'; 
%TARG='AG'; 
%TARG='mPFC';
% ----load channels locations
load('original_chanlocs.mat')
% ----load session info
[~,sessIndx]=xlsread('sessionInfo');

for ii=1:length(subjies)   
    i=subjies{ii};
    da_row=endsWith(sessIndx(:,1), {i});
    % ---load subject data structure
    clearvars dataStruct
    load(['CHODME_' i '.mat']);
    
    for sess=1:3 
        SESSION_NUMBER = sessIndx{da_row,sess+1};
        SESSION_NUMBER = SESSION_NUMBER(2);
              
        % --- Load eeg data (in ft data structure)
        clearvars file_name file_idx eeg_data
        file_name = ['CHODME_' i '_' SESSION_NUMBER '_' TARG...
                 '_TEPs(PRE)_cleaned_ft.mat'];
        file_idx = find(strcmp({dataStruct.filename}, file_name));
        eeg_data = dataStruct(file_idx).data.data_SIR_FT;
        
                if ismember('Afz', eeg_data.label)
                eeg_data.label(strcmpi('Afz', eeg_data.label)) = {'AFz'}; 
                end

                if ismember('Afz', eeg_data.elec.label)
                eeg_data.elec.label(strcmpi('Afz', eeg_data.label)) = {'AFz'};
                end

        trialinfo = eeg_data.trialinfo;        
        cfg = [];
        cfg.reref         = 'yes';
        cfg.refchannel    = {'all'};
        cfg.refmethod     = 'avg';
        cfg.baselinewindow  = [-0.8 -0.01]; 
        EEG_interp = ft_preprocessing(cfg,eeg_data);
        
        cfg=[];
        cfg.trials = find(ismember(trialinfo, 2)); % sham condition
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 45;
        cfg.keeptrials = 'no';
        SOUND_SHAM_PRE = ft_timelockanalysis(cfg, EEG_interp);
        
        cfg.trials = find(ismember(trialinfo, 1)); % active condition
        SOUND_REAL_PRE = ft_timelockanalysis(cfg, EEG_interp);
        
        if sess==1 
            TEPs_SHAM_S1_PRE{ii}=SOUND_SHAM_PRE; 
            TEPs_ACTIVE_S1_PRE{ii}=SOUND_REAL_PRE;
        elseif sess==2 
            TEPs_SHAM_S2_PRE{ii}=SOUND_SHAM_PRE; 
            TEPs_ACTIVE_S2_PRE{ii}=SOUND_REAL_PRE;                            
        elseif sess==3 
            TEPs_SHAM_S3_PRE{ii}=SOUND_SHAM_PRE; 
            TEPs_ACTIVE_S3_PRE{ii}=SOUND_REAL_PRE;         
        end

    end
end

%% SAVE 
% SMA
SMAData.SHAM_S1_PRE = TEPs_SHAM_S1_PRE;
SMAData.ACTIVE_S1_PRE = TEPs_ACTIVE_S1_PRE;
SMAData.SHAM_S2_PRE = TEPs_SHAM_S2_PRE;
SMAData.ACTIVE_S2_PRE = TEPs_ACTIVE_S2_PRE;
SMAData.SHAM_S3_PRE = TEPs_SHAM_S3_PRE;
SMAData.ACTIVE_S3_PRE = TEPs_ACTIVE_S3_PRE;

filepath = '';
filename = 'SMA_TEPs_PRE';
save([filepath filename], 'SMAData', '-v7.3')
%
% AG
AGData.SHAM_S1_PRE = TEPs_SHAM_S1_PRE;
AGData.ACTIVE_S1_PRE = TEPs_ACTIVE_S1_PRE;
AGData.SHAM_S2_PRE = TEPs_SHAM_S2_PRE;
AGData.ACTIVE_S2_PRE = TEPs_ACTIVE_S2_PRE;
AGData.SHAM_S3_PRE = TEPs_SHAM_S3_PRE;
AGData.ACTIVE_S3_PRE = TEPs_ACTIVE_S3_PRE;

filepath = '';
filename = 'AG_TEPs_PRE';
save([filepath filename], 'AGData', '-v7.3')

% mPFC
mPFCData.SHAM_S1_PRE = TEPs_SHAM_S1_PRE;
mPFCData.ACTIVE_S1_PRE = TEPs_ACTIVE_S1_PRE;
mPFCData.SHAM_S2_PRE = TEPs_SHAM_S2_PRE;
mPFCData.ACTIVE_S2_PRE = TEPs_ACTIVE_S2_PRE;
mPFCData.SHAM_S3_PRE = TEPs_SHAM_S3_PRE;
mPFCData.ACTIVE_S3_PRE = TEPs_ACTIVE_S3_PRE;

filepath = '';
filename = 'mPFC_TEPs_PRE';
save([filepath filename], 'mPFCData', '-v7.3')
 %% Grandaverage
% ----Load data
load('mPFC_TEPs_PRE.mat')

% Data = SMAData; 
% Data = AGData;
 Data = mPFCData;

TEPs_SHAM_S1_PRE = Data.SHAM_S1_PRE;
TEPs_ACTIVE_S1_PRE= Data.ACTIVE_S1_PRE;

TEPs_SHAM_S2_PRE = Data.SHAM_S2_PRE;
TEPs_ACTIVE_S2_PRE = Data.ACTIVE_S2_PRE;

TEPs_SHAM_S3_PRE = Data.SHAM_S3_PRE;
TEPs_ACTIVE_S3_PRE = Data.ACTIVE_S3_PRE;

% $$$ SHAM and ACTIVE TEPs 
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
% Session1
GA_TEPs_SHAM_S1_PRE=ft_timelockgrandaverage(cfg,TEPs_SHAM_S1_PRE{:}); 
GA_TEPs_ACTIVE_S1_PRE=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_S1_PRE{:}); 
% Session2
GA_TEPs_SHAM_S2_PRE=ft_timelockgrandaverage(cfg,TEPs_SHAM_S2_PRE{:});
GA_TEPs_ACTIVE_S2_PRE=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_S2_PRE{:});
% Session3
GA_TEPs_SHAM_S3_PRE=ft_timelockgrandaverage(cfg,TEPs_SHAM_S3_PRE{:}); 
GA_TEPs_ACTIVE_S3_PRE=ft_timelockgrandaverage(cfg,TEPs_ACTIVE_S3_PRE{:});

% Average Over Three Sesstions

GA_TEPs_SHAM_ALL_PRE = GA_TEPs_SHAM_S1_PRE;
GA_TEPs_ACTIVE_ALL_PRE = GA_TEPs_ACTIVE_S1_PRE;
for sj = 1:24
   GA_TEPs_SHAM_ALL_PRE.individual(sj,:,:) =  mean( cat(1,GA_TEPs_SHAM_S1_PRE.individual(sj,:,:), GA_TEPs_SHAM_S2_PRE.individual(sj,:,:),GA_TEPs_SHAM_S3_PRE.individual(sj,:,:)), 1);
   GA_TEPs_ACTIVE_ALL_PRE.individual(sj,:,:) = mean( cat(1,GA_TEPs_ACTIVE_S1_PRE.individual(sj,:,:), GA_TEPs_ACTIVE_S2_PRE.individual(sj,:,:), GA_TEPs_ACTIVE_S3_PRE.individual(sj,:,:)) ,1);
end

   
% Cleaned TEPs (Subtraction)
cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
% session1
GA_TEPs_clean_S1_PRE=ft_math(cfg, GA_TEPs_ACTIVE_S1_PRE, GA_TEPs_SHAM_S1_PRE);
% session2
GA_TEPs_clean_S2_PRE=ft_math(cfg, GA_TEPs_ACTIVE_S2_PRE, GA_TEPs_SHAM_S2_PRE);
% session3
GA_TEPs_clean_S3_PRE=ft_math(cfg, GA_TEPs_ACTIVE_S3_PRE, GA_TEPs_SHAM_S3_PRE);
% Average across all sessions
GA_TEPs_clean_ALL_PRE=ft_math(cfg, GA_TEPs_ACTIVE_ALL_PRE, GA_TEPs_SHAM_ALL_PRE);
