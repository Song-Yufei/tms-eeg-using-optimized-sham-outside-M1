% this script calculates Spatial concordance correlation coefficient (CCC) 
% pair-wisely between sessions (S1 v S2, S1 v S3, S2 v S3) and performs
% statistics using one-sample permutation t-tests. 
% require: install matlab function for stat: mult_comp_perm_t1(), author David Groppe
%%  Spatial CCC

% ---- use grand average tep: load_grandAvg_tep.m

% ----load peak latencies info 
load('peak_latency_ranges_sma.mat') 

% #### assign data ####
% sham
COND1 = GA_TEPs_SHAM_S1_PRE;
COND2 = GA_TEPs_SHAM_S2_PRE;
COND3 = GA_TEPs_SHAM_S3_PRE;

% active
% COND1 = GA_TEPs_ACTIVE_S1_PRE;
% COND2 = GA_TEPs_ACTIVE_S2_PRE;
% COND3 = GA_TEPs_ACTIVE_S3_PRE;

% cleaned (active-sham)
% cfg=[];
% cfg.operation =  'subtract';
% cfg.parameter = 'individual';
% COND1 = ft_math(cfg, GA_TEPs_ACTIVE_S1_PRE, GA_TEPs_SHAM_S1_PRE);
% COND2 = ft_math(cfg, GA_TEPs_ACTIVE_S2_PRE, GA_TEPs_SHAM_S2_PRE);
% COND3 = ft_math(cfg, GA_TEPs_ACTIVE_S3_PRE, GA_TEPs_SHAM_S3_PRE);

% #### calculate CCC ####
corrComp = {'c1','c2','c3'};
conditionName = {'S1 & S2', 'S1 & S3', 'S2 & S3'};

time = linspace(-1,0.999,2000)*1000; 

clearvars z ccc avZ rFromZ
% Lin's Concordance Correlation Coefficient
for idx = 1:size(COND1.individual,1)
    for tx = 1:size(COND1.individual,3)    
        ccc.c1(idx,tx) = concordanceCorrelationCoefficient( squeeze(COND1.individual(idx,:,tx))', squeeze(COND2.individual(idx,:,tx))');
        ccc.c2(idx,tx) = concordanceCorrelationCoefficient( squeeze(COND1.individual(idx,:,tx))', squeeze(COND3.individual(idx,:,tx))');
        ccc.c3(idx,tx) = concordanceCorrelationCoefficient( squeeze(COND2.individual(idx,:,tx))', squeeze(COND3.individual(idx,:,tx))');
    end
end

for t = 1:size(ccc.c1,2)
    % Fisher's r to z transform 
    for idx = 1:size(ccc.c1,1)               
        z.c1(idx,t)=.5.*log((1+ccc.c1(idx,t))./(1-ccc.c1(idx,t)));
        z.c2(idx,t)=.5.*log((1+ccc.c2(idx,t))./(1-ccc.c2(idx,t)));
        z.c3(idx,t)=.5.*log((1+ccc.c3(idx,t))./(1-ccc.c3(idx,t)));
    end
      
    % Average of z scores
    avZ.c1(t) = mean(z.c1(:,t),1);
    avZ.c2(t) = mean(z.c2(:,t),1);
    avZ.c3(t) = mean(z.c3(:,t),1);
    
    % Fisher's z to r tranformation for ploting 
    rFromZ.c1(t) = (exp(1)^(2.*avZ.c1(t))-1)/(exp(1)^(2.* avZ.c1(t))+1);
    rFromZ.c2(t) = (exp(1)^(2.*avZ.c2(t))-1)/(exp(1)^(2.* avZ.c2(t))+1);
    rFromZ.c3(t) = (exp(1)^(2.*avZ.c3(t))-1)/(exp(1)^(2.* avZ.c3(t))+1);
            
end

% #### Stat: one sample permutation test
[~,tp1] = min(abs(time - 20));
[~,tp2] = min(abs(time - 300));
ti = tp1:tp2; 

corrP.c1 = mult_comp_perm_t1(z.c1(:,ti),10000,0,0.05); 
corrP.c2 = mult_comp_perm_t1(z.c2(:,ti),10000,0,0.05); 
corrP.c3 = mult_comp_perm_t1(z.c3(:,ti),10000,0,0.05); 

% #### Plot ####
[~,t1] = min(abs(time - -500));
[~,t2] = min(abs(time - -5));
[~,t3] = min(abs(time - 20));
[~,t4] = min(abs(time - 500));

tp = [];
for latx = 1:length(latencies.AVE) % peak time point
    tp(1,latx) = latencies.AVE(latx);
end

fig = figure('color','w');
set(gcf,'position',[250,60,740,580]);
c = get(0, 'DefaultAxesColorOrder');

posA = [0.2,0.77,0.55,0.2];
subplot('position',posA)

for cx = 1:length(conditionName)
    gmfaM  = rFromZ.(corrComp{cx});
    plot(time(t1:t2),gmfaM(t1:t2),'color',c(cx+3,:),'linewidth',2);hold on;
    ps.(['h',num2str(cx)]) = plot(time(t3:t4),gmfaM(t3:t4),'color',c(cx+3,:),'linewidth',2); 
end

% plot peaks
for ix = 1:length(tp)
    plot([tp(ix),tp(ix)],[-5,5],'color',[0.5,0.5,0.5]);
    text(tp(ix),1.1,num2str(tp(ix)),'fontsize',10,'HorizontalAlignment','center');
end

% plot horizantal lines
plot([-500,500],[0,0],'--','color',[0.7,0.7,0.7]);
plot([-500,500],[0.2,0.2],'--','color',[0.7,0.7,0.7]);
plot([-500,500],[0.4,0.4],'--','color',[0.7,0.7,0.7]);
plot([-500,500],[0.6,0.6],'--','color',[0.7,0.7,0.7]); 
plot([-500,500],[0.8,0.8],'--','color',[0.7,0.7,0.7]); 

plot([0,0],[-5,5],'k--','linewidth',1.5);

set(gca,'box','off','xlim',[-50,300],'ylim',[-0.3,1],'tickdir','out','linewidth',1.5,'fontsize',10);
xlabel('Time (ms)');
ylabel('CCC'); % 
yticks([0 0.2 0.4 0.6 0.8 1])

% Plot boxes
rectangle('Position',[20                -1 latRange.AVE(1,2)-20                2],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)
rectangle('Position',[latRange.AVE(2,1) -1 latRange.AVE(2,2)-latRange.AVE(2,1) 2],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)
rectangle('Position',[latRange.AVE(3,1) -1 latRange.AVE(3,2)-latRange.AVE(3,1) 2],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)    
rectangle('Position',[latRange.AVE(4,1) -1 latRange.AVE(4,2)-latRange.AVE(4,1) 2],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)    
rectangle('Position',[latRange.AVE(5,1) -1 300-latRange.AVE(5,1)               2],'LineStyle','--','EdgeColor',[0.7,0.7,0.7],...
    'LineWidth',1)    

% Plot significant time points
timec = time(ti);

logc.c1 = corrP.c1< 0.05;
cc1 = ones(1,length(ti))*-0.15;
cc1(logc.c1==0) = NaN;
plot(timec,cc1,'color',c(4,:),'linewidth',1.2);

logc.c2 = corrP.c2< 0.05;
cc2 = ones(1,length(ti))*-0.2;
cc2(logc.c2==0) = NaN;
plot(timec,cc2,'color',c(5,:),'linewidth',1.2);

logc.c3 = corrP.c3< 0.05;
cc3 = ones(1,length(ti))*-0.25;
cc3(logc.c3==0) = NaN;
plot(timec,cc3,'color',c(6,:),'linewidth',1.2);

lgd = legend([ps.h1,ps.h2,ps.h3],conditionName,'box','off','location','southeast','fontsize',10);
lgd.Position = [0.78,0.92,0.1,0.05];