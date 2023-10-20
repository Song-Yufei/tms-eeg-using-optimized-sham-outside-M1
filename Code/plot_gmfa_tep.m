% This script plot time course of GMFA and topographies of TEPs 
%% 
% ---- load peaks
load('peak_latency_ranges_ag.mat') 

% ---- use grand average gmfa: load_grandAvg_gmfa.m
COND1 = GA_GMFA_SHAM_S1_PRE;
COND2 = GA_GMFA_ACTIVE_S1_PRE;
gmfa.c1 = squeeze(COND1.individual);
gmfa.c2 = squeeze(COND2.individual);

gmfacond = {'c1','c2'};
conditionName = {'Sham','Active'};

time = linspace(-1,0.999,2000)*1000; % follow fieldtrip structure
[~,t1] = min(abs(time - -500));
[~,t2] = min(abs(time - -5));
[~,t3] = min(abs(time - 18));
[~,t4] = min(abs(time - 500));

tp = [];
for latx = 1:length(latencies.AVE) % peak time point
    tp(1,latx) = latencies.AVE(latx);
end

TOI = [];
for latx = 1:size(latRange.AVE,1) % peak time point
    TOI{latx,1}(1,1) = latRange.AVE(latx,1);
    TOI{latx,1}(1,2) = latRange.AVE(latx,2);
end

fig = figure('color','w');
set(gcf,'position',[250,60,740,580]);
c = get(0, 'DefaultAxesColorOrder');

posA = [0.2,0.77,0.55,0.2];
subplot('position',posA)

for cx = 1:length(conditionName)

    gmfaM  = mean(gmfa.(gmfacond{cx}),1);
    gmfaSE = std(gmfa.(gmfacond{cx}),[],1)./sqrt(size(gmfa.(gmfacond{cx}),1));

    plot(time(t1:t2),gmfaM(t1:t2),'color',c(cx,:),'linewidth',2); hold on;
    f = fill([time(t1:t2),fliplr(time(t1:t2))],[gmfaM(t1:t2)-gmfaSE(t1:t2),fliplr(gmfaM(t1:t2)+gmfaSE(t1:t2))],c(cx,:));
    set(f,'FaceAlpha',0.2);set(f,'EdgeColor', 'none');

    ps.(['h',num2str(cx)]) = plot(time(t3:t4),gmfaM(t3:t4),'color',c(cx,:),'linewidth',2); hold on;
    f = fill([time(t3:t4),fliplr(time(t3:t4))],[gmfaM(t3:t4)-gmfaSE(t3:t4),fliplr(gmfaM(t3:t4)+gmfaSE(t3:t4))],c(cx,:));
    set(f,'FaceAlpha',0.2);set(f,'EdgeColor', 'none');
end

% plot peaks
for ix = 1:length(tp)
    plot([tp(ix),tp(ix)],[-5,6],'color',[0.5,0.5,0.5]);
    text(tp(ix),6.4,num2str(tp(ix)),'fontsize',10,'HorizontalAlignment','center');

end

plot([0,0],[-5,6],'k--','linewidth',1.5);

set(gca,'box','off','xlim',[-50,300],'ylim',[0,6],'tickdir','out','linewidth',1.5,'fontsize',10);
xlabel('Time (ms)');
ylabel('GMFA (\muV)'); 
% Plot boxes
rectangle('Position',[18                -1 latRange.AVE(1,2)-18                7],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)
rectangle('Position',[latRange.AVE(2,1) -1 latRange.AVE(2,2)-latRange.AVE(2,1) 7],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)
rectangle('Position',[latRange.AVE(3,1) -1 latRange.AVE(3,2)-latRange.AVE(3,1) 7],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)    
rectangle('Position',[latRange.AVE(4,1) -1 latRange.AVE(4,2)-latRange.AVE(4,1) 7],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)    
rectangle('Position',[latRange.AVE(5,1) -1 300-latRange.AVE(5,1)               7],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)     

lgd = legend([ps.h1,ps.h2],conditionName,'box','off','location','southeast','fontsize',10);
lgd.Position = [0.76,0.92,0.1,0.05];
%% 
% ---- load peaks and channel locations
load('original_chanlocs.mat');
load('peak_latency_ranges_ag.mat') 

% ---- use grand average tep: load_grandAvg_tep.m
COND1 = GA_TEPs_SHAM_ALL_PRE;
COND2 = GA_TEPs_ACTIVE_ALL_PRE;

cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'individual';
COND3 =ft_math(cfg, COND2, COND1);

time = linspace(-1,0.999,2000)*1000; 

% peaks
tp = [];
for latx = 1:length(latencies.AVE) 
    tp(1,latx) = latencies.AVE(latx);
end
 
 %#### Topoplots ####
fig = figure('color','w');
set(gcf,'position',[250,60,740,580]);

twidth = 0.11;
theight = 0.11;

txpos = linspace(0.15,0.7,length(tp));
typos = [0.56,0.4,0.24]; 
ix = 4:length(tp)*3+3;
ix = reshape(ix,3,[]);

for plotx = 1:length(tp)
    
    [~,tInd2plot] = min(abs(time-tp(plotx)));
    
    % CON1
    posName = ['pos',num2str(ix(1,plotx))];
    pos.(posName) = [txpos(plotx),typos(1),twidth,theight];
    subplot('position',pos.(posName))
    topoplotIndie(mean(COND1.individual(:,:,tInd2plot),1),original_chanlocs,'electrodes','off','numcontour',0);
    set(gca,'clim',[-1 1]*2)
    text(0,0.7,[num2str(tp(plotx)),'ms'],'horizontalalignment','center','fontsize',10);

    if plotx == length(tp)
        c = colorbar;
        tmpc = c.Position;
        c.Position = [tmpc(1)+0.07,tmpc(2),tmpc(3)+0.005,0.05];
        c.LineWidth = 1.5;
        c.Ticks = [-2;0;2];
        title(c,'\muV');
    end  

    % CON2
    posName = ['pos',num2str(ix(2,plotx))];
    pos.(posName) = [txpos(plotx),typos(2),twidth,theight];
    subplot('position',pos.(posName))
    topoplotIndie(mean(COND2.individual(:,:,tInd2plot),1),original_chanlocs,'electrodes','off','numcontour',0);
    set(gca,'clim',[-1 1]*2)

    if plotx == length(tp)
        c = colorbar;
        tmpc = c.Position;
        c.Position = [tmpc(1)+0.07,tmpc(2),tmpc(3)+0.005,0.05];
        c.LineWidth = 1.5;
        c.Ticks = [-2;0;2];
        title(c,'\muV');
    end  
    
    % COND 3
    posName = ['pos',num2str(ix(3,plotx))];
    pos.(posName) = [txpos(plotx),typos(3),twidth,theight];
    subplot('position',pos.(posName))
    topoplotIndie(mean(COND3.individual(:,:,tInd2plot),1),original_chanlocs,'electrodes','off','numcontour',0);
    set(gca,'clim',[-1 1]*2)

    if plotx == length(tp)
        c = colorbar;
        tmpc = c.Position;
        c.Position = [tmpc(1)+0.07,tmpc(2),tmpc(3)+0.005,0.05];
        c.LineWidth = 1.5;
        c.Ticks = [-2;0;2];
        title(c,'\muV');
    end  
end

[map, descriptorname, description] = colorcet('Gouldian');
colormap(map)
