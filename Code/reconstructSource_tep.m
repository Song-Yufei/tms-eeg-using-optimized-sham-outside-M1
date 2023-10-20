% This script loads the preprocessed dataset after the tms-eeg cleaning pipeline, 
% and reconstructs source activity of TEPs using MNE method.
% Note: since leadfield matrix was created without bad channels. 
% The interpolated channels have to be removed for inverse estimation.
%% load data and MNE inverse estimation

subjies={'002','003', '004', '007', '009', '010', '011', '012',...
    '013', '014', '015', '017', '018', '019', '020', '021', '022','024','025','027','028','030','031','032'};

TARG={'mPFC','AG','SMA'};
% ----load session info
[~,sessIndx]=xlsread('sessionInfo');

for Ti = 1:length(TARG)
 
    for ii=1:length(subjies)  
    i=subjies{ii};
    da_row=endsWith(sessIndx(:,1), {i});
    % ---load subject data structure
    clearvars dataStruct
    load(['CHODME_' i '.mat']);
    
            for sess=1:3 
                
                SESSION_NUMBER = sessIndx{da_row,sess+1};
                SESSION_NUMBER = SESSION_NUMBER(2);
                 
                clearvars Mdw indsp LF LnrN headmodel... 
                 source_TEPs_REAL_PRE source_TEPs_SHAM_PRE...
                 eeg_data lf_label 

                 % ---load leadfield matrix
                 clearvars file_name file_idx 
                 file_name = ['CHODME_' i '_' SESSION_NUMBER '_' TARG{Ti}...
                     '_leadfield.mat'];
                 file_idx = find(strcmp({dataStruct.filename}, file_name));
                 lf_label = dataStruct(file_idx).data.LF.label;
                 LnrN = dataStruct(file_idx).data.LnrN;
                 
                 % --- load individual headmodel
                 clearvars file_name file_idx 
                 file_name = 'headmodel.mat';
                 file_idx = find(strcmp({dataStruct.filename}, file_name));
                 headmodel = dataStruct(file_idx).data.headmodel;           
                 
                 % --- Load eeg data (in ft data structure)
                 clearvars file_name file_idx 
                 file_name = ['CHODME_' i '_' SESSION_NUMBER '_' TARG{Ti}...
                     '_TEPs(PRE)_cleaned_ft.mat'];
                 file_idx = find(strcmp({dataStruct.filename}, file_name));
                 eeg_data = dataStruct(file_idx).data.data_SIR_FT;
                
                % 
                if ismember('Afz', eeg_data.label)
                eeg_data.label(strcmpi('Afz', eeg_data.label)) = {'AFz'}; 
                end

                if ismember('Afz', eeg_data.elec.label)
                eeg_data.elec.label(strcmpi('Afz', eeg_data.label)) = {'AFz'};
                end

%#### compute inverse estimation with MNE ####
% calulate source model matrix using individual headmodel
distanceTH=5; 
attenuation_length=10; 
[Mdw,indsp]=createSourceCovM(headmodel, distanceTH, attenuation_length);

%average trials of interest 
clearvars eeg_data_SOURCE eeg_data_AVG Xst;
cfg = [];
cfg.keeptrials = 'no';
cfg.channel = lf_label';% remove channels interpolated
cfg.trials=find(eeg_data.trialinfo(:,1)==1 ); % Active condition
eeg_data_AVG= ft_timelockanalysis(cfg, eeg_data);
Xst=eeg_data_AVG.avg;
[eeg_data_SOURCE, ~, ~]=estimateCorrelatedSourceAmplitudesFast(LnrN, Mdw, Xst, 0.01, indsp, 1:size(Xst,2), [1 size(Xst,2)]);
source_TEPs_REAL_PRE=eeg_data_SOURCE;

clearvars eeg_data_SOURCE eeg_data_AVG Xst;
cfg = [];
cfg.keeptrials = 'no';
cfg.channel = lf_label';
cfg.trials=find(eeg_data.trialinfo(:,1)==2 ); % Sham condition
eeg_data_AVG= ft_timelockanalysis(cfg, eeg_data);
Xst=eeg_data_AVG.avg;
[eeg_data_SOURCE, ~, ~]=estimateCorrelatedSourceAmplitudesFast(LnrN, Mdw, Xst, 0.01, indsp, 1:size(Xst,2), [1 size(Xst,2)]);
source_TEPs_SHAM_PRE=eeg_data_SOURCE;

%####   
                if sess==1 
                    source_SHAM_S1_PRE{ii}=source_TEPs_SHAM_PRE; 
                    source_REAL_S1_PRE{ii}=source_TEPs_REAL_PRE; 
                elseif sess==2 
                    source_SHAM_S2_PRE{ii}=source_TEPs_SHAM_PRE;  
                    source_REAL_S2_PRE{ii}=source_TEPs_REAL_PRE; 
                elseif sess==3 
                    source_SHAM_S3_PRE{ii}=source_TEPs_SHAM_PRE;
                    source_REAL_S3_PRE{ii}=source_TEPs_REAL_PRE; 
                end

            end 
    end 
            
        if Ti == 1 % mPFC          
            mPFCData_TEPs.SHAM_S1_PRE = source_SHAM_S1_PRE;
            mPFCData_TEPs.ACTIVE_S1_PRE = source_REAL_S1_PRE;
            mPFCData_TEPs.SHAM_S2_PRE = source_SHAM_S2_PRE;
            mPFCData_TEPs.ACTIVE_S2_PRE = source_REAL_S2_PRE;
            mPFCData_TEPs.SHAM_S3_PRE = source_SHAM_S3_PRE;
            mPFCData_TEPs.ACTIVE_S3_PRE = source_REAL_S3_PRE;
            filepath = '';
            filename = 'mPFC_TEPs_Source_PRE';
            save([filepath filename], 'mPFCData_TEPs', '-v7.3')
            
        elseif Ti == 2 % AG               
            AGData_TEPs.SHAM_S1_PRE = source_SHAM_S1_PRE;
            AGData_TEPs.ACTIVE_S1_PRE = source_REAL_S1_PRE;
            AGData_TEPs.SHAM_S2_PRE = source_SHAM_S2_PRE;
            AGData_TEPs.ACTIVE_S2_PRE = source_REAL_S2_PRE;
            AGData_TEPs.SHAM_S3_PRE = source_SHAM_S3_PRE;
            AGData_TEPs.ACTIVE_S3_PRE = source_REAL_S3_PRE;
            filepath = '';
            filename = 'AG_TEPs_Source_PRE';
            save([filepath filename], 'AGData_TEPs', '-v7.3') 
         
          elseif  Ti == 3 %SMA                
            SMAData_TEPs.SHAM_S1_PRE = source_SHAM_S1_PRE;
            SMAData_TEPs.ACTIVE_S1_PRE = source_REAL_S1_PRE;
            SMAData_TEPs.SHAM_S2_PRE = source_SHAM_S2_PRE;
            SMAData_TEPs.ACTIVE_S2_PRE = source_REAL_S2_PRE;
            SMAData_TEPs.SHAM_S3_PRE = source_SHAM_S3_PRE;
            SMAData_TEPs.ACTIVE_S3_PRE = source_REAL_S3_PRE;
            filepath = '';
            filename = 'SMA_TEPs_Source_PRE';
            save([filepath filename], 'SMAData_TEPs', '-v7.3')
        
        end              
        clearvars source_SHAM_S1_PRE source_REAL_S1_PRE source_SHAM_S2_PRE source_REAL_S2_PRE...
            source_SHAM_S3_PRE source_REAL_S3_PRE
end % 

%% Normalize using Z-score Transformation 
% the source activity is standardized using a z-score transformation 
% with respect to Pre-stimulus time (-600 to-100ms) 
% Source map standardization does not alter the dynamics of the source time series and 
% only scales their respective amplitude changes

% Load data
load('AG_TEPs_Source_PRE.mat') % AG, SMA, mPFC

% Data = SMAData_TEPs; 
 Data = AGData_TEPs;
% Data = mPFCData_TEPs;

source_SHAM_S1_PRE = Data.SHAM_S1_PRE;
source_ACTIVE_S1_PRE= Data.ACTIVE_S1_PRE;
source_SHAM_S2_PRE = Data.SHAM_S2_PRE;
source_ACTIVE_S2_PRE = Data.ACTIVE_S2_PRE;
source_SHAM_S3_PRE = Data.SHAM_S3_PRE;
source_ACTIVE_S3_PRE = Data.ACTIVE_S3_PRE;

source_SHAM_PRE = NaN(24,size(source_SHAM_S1_PRE{1,1},1),size(source_SHAM_S1_PRE{1,1},2) );
source_ACTIVE_PRE = NaN(24,size(source_ACTIVE_S1_PRE{1,1},1),size(source_ACTIVE_S1_PRE{1,1},2) );

for sj =1:24
    clearvars ss dummy
    % average across 3 sessions and  then Normalize (sham) 
    ss(1,:,:) = source_SHAM_S1_PRE{1,sj};  
    ss(2,:,:) = source_SHAM_S2_PRE{1,sj};
    ss(3,:,:) = source_SHAM_S3_PRE{1,sj};
    dummy = squeeze(nanmean(ss,1));
    source_SHAM_PRE(sj,:,:)  = (dummy - nanmean(dummy(:,400:900),2))./nanstd(dummy(:,400:900),[],2); 

    clearvars sa dummy
     % average across 3 sessions and  then Normalize (active) 
    sa(1,:,:) = source_ACTIVE_S1_PRE{1,sj};  
    sa(2,:,:) = source_ACTIVE_S2_PRE{1,sj};
    sa(3,:,:) = source_ACTIVE_S3_PRE{1,sj};
    dummy = squeeze(nanmean(sa,1));
    source_ACTIVE_PRE(sj,:,:) = (dummy - nanmean(dummy(:,400:900),2))./nanstd(dummy(:,400:900),[],2);
end

%% Plot Cleaned (Active-Sham) TEPs 
% load peaks
load('peak_latency_ranges_ag.mat')
% load relevant data
load('generic_headmodel.mat')

% load ROI, pre-defined 
load('poi_mpfc.mat') 
load('poi_ag.mat')
load('poi_sma.mat')

% time window
TOI = [];
for latx = 1:size(latRange.AVE,1) % peak time point
    TOI{latx,1}(1,1) = latRange.AVE(latx,1);
    TOI{latx,1}(1,2) = latRange.AVE(latx,2);
end

poi.indx1 = poi_mpfc; % mpfc
poi.indx2 = poi_ag; % left ag
poi.indx3 = poi_sma; % sma

clearvars cond1 cond2 cond3 COND1 COND2 COND3

COND1 = squeeze(mean(source_SHAM_PRE,1));
COND2 = squeeze(mean(source_ACTIVE_PRE,1));
COND3 = COND2 - COND1; % Cleaned (Active-Sham) TEPs 

% time boundary
time =   linspace(-1,0.9990,2000)*1000;
times = time(1001:3:1401);
[~,t1] = min(abs(time - -100));
[~,t2] = min(abs(time - -5));
[~,t3] = min(abs(time - 20));
[~,t4] = min(abs(time - 500));

% color
[map, descriptorname, description] = colorcet('BWRA'); 
z_lim=[-10 10]; 

% 
fig = figure('color','w');
set(gcf,'position',[0.1 0.8 2 0.4]*1000);
c = get(0, 'DefaultAxesColorOrder');
f1 = tiledlayout(2,8,'TileSpacing','compact');

for i =1:5
su(i) = tiledlayout(f1,1,1,'TileSpacing','none');
su(i).Layout.Tile = i+3; 
su(i).Layout.TileSpan = [1 1];
end

for i =1:5
sd(i) = tiledlayout(f1,1,1,'TileSpacing','none');
sd(i).Layout.Tile = i+11; 
sd(i).Layout.TileSpan = [1 1];
end

for i =1:3
td(i) = tiledlayout(f1,1,1,'TileSpacing','none');
td(i).Layout.Tile = i; 
td(i).Layout.TileSpan = [1 1];
end

% #### temporal plot ####
for itp =1:3
    clear poi_indx cond gmfaM gmfaSE
    if itp ==1
    poi_indx = poi.indx1; % mpfc
    elseif itp ==2
    poi_indx = poi.indx3; % sma
    else
    poi_indx = poi.indx2; % ag (left) 
    end 

    cond = squeeze(nanmean(source_ACTIVE_PRE(:,poi_indx,:),2)) - squeeze(nanmean(source_SHAM_PRE(:,poi_indx,:),2));     
    gmfaM = mean(cond,1);
    gmfaSE = std(cond,[],1)./sqrt(size(cond,1));

    nexttile(td(itp))
    plot(time(t1:t2),gmfaM(t1:t2),'color',c(itp,:),'linewidth',2); hold on;
    f = fill([time(t1:t2),fliplr(time(t1:t2))],[gmfaM(t1:t2)-gmfaSE(t1:t2),fliplr(gmfaM(t1:t2)+gmfaSE(t1:t2))],c(itp,:));
    set(f,'FaceAlpha',0.3);set(f,'EdgeColor', 'none');
    
    pg.(['h',num2str(2)]) = plot(time(t3:t4),gmfaM(t3:t4),'color',c(itp,:),'linewidth',2); hold on;
    f = fill([time(t3:t4),fliplr(time(t3:t4))],[gmfaM(t3:t4)-gmfaSE(t3:t4),fliplr(gmfaM(t3:t4)+gmfaSE(t3:t4))],c(itp,:));
    set(f,'FaceAlpha',0.3);set(f,'EdgeColor', 'none');

    hold on,
    plot([0,0],[-10,10],'k--','linewidth',1);
    plot([-100,500],[0,0],'k--','linewidth',0.5);
    set(gca,'box','off','xlim',[-100,500],'ylim',[-10,10],'tickdir','out','linewidth',1,'fontsize',10);
    xlabel('Time (ms)');
    ylabel('z value');
    title('Real-Sham', 'Units', 'normalized', 'Position', [0.8, 0.9, 1],'fontsize',10);
end 
 
%#### spatial plot #####
% average within each TOI
for b=1:size(TOI,1)
   
    [~,temtx1] = min(abs(time-round(TOI{b}(1))));
    [~,temtx2] = min(abs(time-round(TOI{b}(2))));    
    nexttile(su(b));   
    PlotDataOnMesh(headmodel.smesh,mean(COND3(:,temtx1:temtx2),2),'caxis', z_lim,'colormap', map, 'colorbar', 0, 'view', [-110 40]) % view, [-125 35]  [-5 90]
end

for b=1:size(TOI,1)
   
    [~,temtx1] = min(abs(time-round(TOI{b}(1))));
    [~,temtx2] = min(abs(time-round(TOI{b}(2))));    
    nexttile(sd(b));
    PlotDataOnMesh(headmodel.smesh,mean(COND3(:,temtx1:temtx2),2),'caxis', z_lim,'colormap', map, 'colorbar', 0, 'view', [-5 90]) % view, [-125 35] [-110 40]
end

cb = colorbar;
cb.Layout.Tile = "east";
cb.LineWidth = 1;
title(cb,'z','fontsize',12);
