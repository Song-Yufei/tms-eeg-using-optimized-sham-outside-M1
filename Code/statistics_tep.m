% This script compares sham vs. active for each target
% using cluster-based permutation statistics as implemented in FieldTrip. 
%%
% ---- load peaks and neighbours
load('neighbour_layout64_custom.mat')
load('peak_latency_ranges_sma.mat') 

% ---- use grand average tep: load_grandAvg_tep.m
COND1 = GA_TEPs_SHAM_ALL_PRE;
COND2 = GA_TEPs_ACTIVE_ALL_PRE;

TOI = [];
for latx = 1:size(latRange.AVE,1) % peak time point
    TOI{latx,1}(1,1) = latRange.AVE(latx,1)/1000;
    TOI{latx,1}(1,2) = latRange.AVE(latx,2)/1000;
end

subj=size(COND1.individual,1);

% go through each TOI
for b = 1:size(TOI,1)
    cfg = [];
    cfg.latency = TOI{b}; 
    cfg.avgovertime = 'yes';
    cfg.method = 'montecarlo'; 
    cfg.statistic = 'ft_statfun_depsamplesT';
    cfg.correctm = 'cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2; 
    cfg.neighbours = neighbour_layout.custneighbor.neighbours; 
    cfg.tail = 0; 
    cfg.clustertail = 0;
    cfg.correcttail = 'alpha';  
    cfg.alpha = 0.05;
    cfg.numrandomization =2000;
    design = zeros(2,2*subj);
    for i = 1:subj
      design(1,i) = i;
    end
    for i = 1:subj
      design(1,subj+i) = i;
    end
    design(2,1:subj)        = 1;
    design(2,subj+1:2*subj) = 2;
    cfg.design = design;
    cfg.uvar  = 1;
    cfg.ivar  = 2;
    stat = ft_timelockstatistics(cfg, COND2, COND1);    
    statOut{b} = stat;
end
 
save ([],'statOut');
%% PLOT spatial distribution of T-values over the significant cluster
% ---- load statOut
b= 4; % TOI index
stat = statOut{1,b};

pval_lim=0.005; % two-tailed (0.05/2 = 0.025), Bonferroni correction by the number of TOIs (0.025/5=0.005)

fig = figure('color','w');
    cfg = [];
    cfg.comment = [num2str(TOI{b}(1)) '-' num2str(TOI{b}(2)-1) 'ms'];
    cfg.commentpos = 'lefttop';
    cfg.style = 'straight';
    cfg.alpha  = pval_lim;
    cfg.parameter = 'stat';
    cfg.zlim   = [-4 4]; 
    cfg.highlightsymbolseries  = ['.','.','.','.','.'];
    cfg.highlightsizeseries =[2 2 2 ].*10;
    cfg.highlightcolorneg=[	221,125,30]./255;
    cfg.highlightcolorpos=[	221,125,30]./255;
    cfg.subplotsize = [1 5];
    cfg.layout = neighbour_layout.layout;
    [mymap, descriptorname, description] = colorcet('Gouldian');
    cfg.colormap = mymap;
    ft_clusterplot(cfg, stat);        
%% Plot time course 
% ---- load statOut

b= 1;  % indx of TOIs
stat = statOut{1,b};

% ----postive cluster
%chan=find(stat.posclusterslabelmat ==1);

% ----negtive cluster
chan=find(stat.negclusterslabelmat ==1);

% Condition name
condName = {'SHAM','ACTIVE'};

[~,t1] = min(abs(COND1.time*1000 - -500));
[~,t2] = min(abs(COND1.time*1000 - -5));
[~,t3] = min(abs(COND1.time*1000 - 18));
[~,t4] = min(abs(COND1.time*1000 - 500));

time = COND1.time*1000;

c = get(0, 'DefaultAxesColorOrder');

for cx = 1:2
    if cx ==1
    COND = squeeze(mean(COND1.individual(:,chan,:),2)); %sham
    else
    COND = squeeze(mean(COND2.individual(:,chan,:),2)); %active
    end
    
    gmfaM = mean(COND,1);
    gmfaSE = std(COND,[],1)./sqrt(size(COND,1));
    
    % %positive
    %subplot('position',[0.45,0.6,0.3,0.25])
    % negtive
    subplot('position',[0.45,0.25,0.3,0.25])

    plot(time(t1:t2),gmfaM(t1:t2),'color',c(cx,:),'linewidth',2); hold on;
    f = fill([time(t1:t2),fliplr(time(t1:t2))],[gmfaM(t1:t2)-gmfaSE(t1:t2),fliplr(gmfaM(t1:t2)+gmfaSE(t1:t2))],c(cx,:));
    set(f,'FaceAlpha',0.3);set(f,'EdgeColor', 'none');
    
    pg.(['h',num2str(cx)]) = plot(time(t3:t4),gmfaM(t3:t4),'color',c(cx,:),'linewidth',2); hold on;
    f = fill([time(t3:t4),fliplr(time(t3:t4))],[gmfaM(t3:t4)-gmfaSE(t3:t4),fliplr(gmfaM(t3:t4)+gmfaSE(t3:t4))],c(cx,:));
    set(f,'FaceAlpha',0.3);set(f,'EdgeColor', 'none');
   clear COND 
end

plot([0,0],[-5,5],'k--','linewidth',2);

set(gca,'box','off','xlim',[-50,300],'ylim',[-5,5],'tickdir','out','linewidth',2,'fontsize',12);
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');

lgd1 = legend([pg.h1,pg.h2],condName,'box','off','location','southeast','fontsize',10);
%lgd1.Position = [0.66,0.63,0.1,0.05]; %positive
lgd1.Position = [0.66,0.28,0.1,0.05]; %negative
