% this script calculates temporal concordance correlation coefficient (CCC) 
% pair-wisely between sessions (S1 v S2, S1 v S3, S2 v S3) and performs
% statistics using one-sample permutation t-tests. 
% install matlab function for stat: mult_comp_perm_t1(), author David Groppe
%%  Temporal CCC 

% ---- use grand average tep: load_grandAvg_tep.m

% ----load neighbours for plot
load('neighbour_layout64_custom.mat')

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
time = linspace(-1,0.999,2000)*1000;% in ms
corrComp = {'c1','c2','c3'};

% ----load 5 TOIs
 clear TOI
 load('peak_latency_ranges_sma.mat') 
for latx = 1:size(latRange.AVE,1) 
    TOI{latx,1}(1,1) = latRange.AVE(latx,1);
    TOI{latx,1}(1,2) = latRange.AVE(latx,2);
end     

% calculate ccc (Lin's Concordance Correlation Coefficient)
clearvars z ccc avZ rFromZ
for idx = 1:size(COND1.individual,1)
    for chx = 1:size(COND1.individual,2)    

        for b = 1:size(TOI,1)+1
            if b<6 
            [~,temtx1] = min(abs(time-round(TOI{b}(1))));
            [~,temtx2] = min(abs(time-round(TOI{b}(2))));             
            
            ccc.c1{b}(idx,chx) = concordanceCorrelationCoefficient( squeeze(COND1.individual(idx,chx,temtx1:temtx2-1)), squeeze(COND2.individual(idx,chx,temtx1:temtx2-1)) );
            ccc.c2{b}(idx,chx) = concordanceCorrelationCoefficient( squeeze(COND1.individual(idx,chx,temtx1:temtx2-1)), squeeze(COND3.individual(idx,chx,temtx1:temtx2-1)) );
            ccc.c3{b}(idx,chx) = concordanceCorrelationCoefficient( squeeze(COND2.individual(idx,chx,temtx1:temtx2-1)), squeeze(COND3.individual(idx,chx,temtx1:temtx2-1)) ); 
           
            else  
            [~,temtx1] = min(abs(time-20));
            [~,temtx2] = min(abs(time-300));    
            ccc.c1{b}(idx,chx) = concordanceCorrelationCoefficient( squeeze(COND1.individual(idx,chx,temtx1:temtx2)), squeeze(COND2.individual(idx,chx,temtx1:temtx2)) );
            ccc.c2{b}(idx,chx) = concordanceCorrelationCoefficient( squeeze(COND1.individual(idx,chx,temtx1:temtx2)), squeeze(COND3.individual(idx,chx,temtx1:temtx2)) );
            ccc.c3{b}(idx,chx) = concordanceCorrelationCoefficient( squeeze(COND2.individual(idx,chx,temtx1:temtx2)), squeeze(COND3.individual(idx,chx,temtx1:temtx2)) ); 
  
           end
        end
    end
end

% Fisher's r to z transformation
for ch = 1:size(ccc.c1{1},2)
    for idx = 1:size(ccc.c1{1},1)
       for b = 1:size(TOI,1)+1  
        z.c1{b}(idx,ch)=.5.*log((1+ccc.c1{b}(idx,ch))./(1-ccc.c1{b}(idx,ch)));
        z.c2{b}(idx,ch)=.5.*log((1+ccc.c2{b}(idx,ch))./(1-ccc.c2{b}(idx,ch)));
        z.c3{b}(idx,ch)=.5.*log((1+ccc.c3{b}(idx,ch))./(1-ccc.c3{b}(idx,ch)));      
       end
    end
end    


for ch = 1:size(ccc.c1{1},2)
    
    for b =1:size(TOI,1)+1
        
        %  Average of z scores
        avZ.c1{b}(ch) = mean(z.c1{b}(:,ch),1);
        avZ.c2{b}(ch) = mean(z.c2{b}(:,ch),1);
        avZ.c3{b}(ch) = mean(z.c3{b}(:,ch),1);
        
        % Transform meansubjects' z back to r for plotting 
        rFromZ.c1{b}(ch) = (exp(1).^(2.*avZ.c1{b}(ch))-1)./(exp(1).^(2.*avZ.c1{b}(ch))+1);
        rFromZ.c2{b}(ch) = (exp(1).^(2.*avZ.c2{b}(ch))-1)./(exp(1).^(2.*avZ.c2{b}(ch))+1);
        rFromZ.c3{b}(ch) = (exp(1).^(2.*avZ.c3{b}(ch))-1)./(exp(1).^(2.*avZ.c3{b}(ch))+1);
                
    end
    
end

% #### Stat: one sample permutation test
for b =1:size(TOI,1)+1   
    corrP.c1{b} = mult_comp_perm_t1(z.c1{b},10000, 0, 0.05); 
    corrP.c2{b} = mult_comp_perm_t1(z.c2{b},10000, 0, 0.05); 
    corrP.c3{b} = mult_comp_perm_t1(z.c3{b},10000, 0, 0.05); 
end

for b =1:size(TOI,1)+1
    logcc.c1{b} = corrP.c1{b}<0.05; 
    logcc.c2{b} = corrP.c2{b}<0.05;
    logcc.c3{b} = corrP.c3{b}<0.05;
end 

% #### PLot CCC with significant channels####
fig = figure('color','w');
set(gcf,'position',[250,60,740,580]);

for cx = 1: length(corrComp)
   for b = 1:6    
        HighChan = [];
        plotStruc = [];
        plotStruc.time = 1; 
        plotStruc.dimord = 'chan';
        plotStruc.label = COND1.label;
        plotStruc.avg = rFromZ.(corrComp{cx}){b};  
        HighChan = logcc.(corrComp{cx}){b}(1,:); 

        if cx == 1
            locat = b+6;
        elseif cx ==2
            locat = b+12;
        else
            locat = b+18;
        end

        cfg = [];
        cfg.figure     = subplot(4,6,locat);
        cfg.layout = neighbour_layout.layout;
        [mymap, descriptorname, description] = colorcet('D1A');
        cfg.colormap = mymap;
        cfg.interactive = 'no';
        cfg.zlim = [-1,1];
        cfg.comment = 'no';
        cfg.style= 'straight';
        cfg.marker = 'off';
        cfg.highlight = 'on';
        cfg.highlightchannel   =  COND1.label(HighChan); 
        cfg.highlightsymbol    = 'x';
        cfg.highlightsize = 4;
        cfg.highlightcolor   = [0 0 0];  
        ft_topoplotER(cfg,plotStruc);
        
        if locat ==24
        c = colorbar;
        tmpc = c.Position;
        c.Position = [tmpc(1)+0.07,tmpc(2)-0.028,tmpc(3)+0.005,0.55];
        c.LineWidth = 1.5;
        c.Ticks = [-1;-0.8;-0.6;-0.4;-0.2;0;0.2;0.4;0.6;0.8;1];
        c.TickDirection = 'out';
        title(c,'CCC');
        end            
   end
end