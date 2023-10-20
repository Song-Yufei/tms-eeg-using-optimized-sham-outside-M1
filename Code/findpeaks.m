% Find peaks
% This script is adapted from the code:
% https://github.com/nigelrogasch/DXM_TMS-EEG_paper/tree/master/find_peaks
% Author: Nigel Rogasch, Monash University
% It uses averaged GA_GMFAs across the sham and active conditions and automatically 
% finds peaks per target. For target AG and mPFC in which the
% early peaks are not identifiable, the peaks are determined using evoked potentials,
% averaged across electrodes near the target.
% peak latencies and ranges used in the current paper are saved in the supportFunctionAndFile
% folder.
%%
% ---- use grand average gmfa: load_grandAvg_gmfa.m

% two different dataset (REAL is computed for visulization purpose)
cond = {'REAL';'AVE'}; 

% Time window to search for peaks (18-300 ms)
peakWin = [18,300];

% Number of +/- ms for peak
peakDef = [5,15];
peakDefWin = [13,100;101,300]; 
time = GA_GMFA_ACTIVE_ALL_PRE.time*1000;
    
gmfaAve.(cond{1,1}) = squeeze(mean(GA_GMFA_ACTIVE_ALL_PRE.individual,1))';
gmfaAve.(cond{2,1}) = squeeze(mean(GA_GMFA_SHAM_ALL_PRE.individual,1))';
%(sham+active)/2
gmfaAve.(cond{2,1}) = mean( cat(1,gmfaAve.(cond{1,1}),gmfaAve.(cond{2,1})),1 ); 

% Automatically find GMFA peaks
% Peak is defined as a point which is larger in amplitude that the point
% +/- 5 ms.
for sidx = 1:size(cond,1)
    % Find time points to define time series
    [val,tpW(1,1)] = min(abs(time-peakWin(1,1)));
    [val,tpW(1,2)] = min(abs(time-peakWin(1,2)));
    % Extract time series
    tseries = gmfaAve.(cond{sidx,1});
    % Set -5 to 18 ms to 0;
    [~,tp1] = min(abs(-5-time));
    [~,tp2] = min(abs(17-time));
    tseries(1,tp1:tp2) = 0;
    latHold = [];
    num = 1;
    for b = tpW(1,1):tpW(1,2)    
        ztime = time(1,b);        
        % Find +/- windows
        tPlus = [];
        tMinus = [];
        if time(1,b) <= peakDefWin(1,2)
            for c = 1:peakDef(1,1)
                tPlus(c,1) = tseries(1,b) - tseries(1,b+c);
                tMinus(c,1) = tseries(1,b) - tseries(1,b-c);
            end
        else
            for c = 1:peakDef(1,2)
                tPlus(c,1) = tseries(1,b) - tseries(1,b+c);
                tMinus(c,1) = tseries(1,b) - tseries(1,b-c);
            end
        end

        % Find time points greater than 0
        tPlusLog = tPlus > 0;
        tMinusLog = tMinus > 0;        
        testOut(1,b-(tpW(1,1)-1)) = size(tPlus,1) + size(tMinus,1);

        if  time(1,b) <= peakDefWin(1,2)
            if sum(tPlusLog) + sum(tMinusLog) == peakDef(1,1)*2
                latHold(num,1) = b;
                num = num+1;
            end
        else
            if sum(tPlusLog) + sum(tMinusLog) == peakDef(1,2)*2
                latHold(num,1) = b;
                num = num+1;
            end
        end
    end
    
    % Calculate latencies
    latencies.(cond{sidx,1}) = time(1,latHold);
end

%#### Find missing early peaks using grand averaged evoked potentials #####
% take (sham+active)/2
%{
for sj = 1:24
baseMean(sj,:,:) = mean(cat(1, GA_TEPs_SHAM_ALL_PRE.individual(sj,:,:),GA_TEPs_ACTIVE_ALL_PRE.individual(sj,:,:)), 1);
end
dataMean = squeeze(mean(baseMean,1));
%}

peakWinOrig = peakWin;
condOrig = cond;
gmfaOrig = gmfaAve;

peakWin = [20,90];
cond = {'AVE'};

clearvars chan latenciesTemp gmfaAve.AVE
% chan=ismember(GA_TEPs_ACTIVE_ALL_PRE.label, {'FCZ' 'Cz' 'Fz' 'FC1' 'FC2'}); % sma
% chan=ismember(GA_TEPs_ACTIVE_ALL_PRE.label, {'AFz' 'AF3' 'AF4' 'F1' 'F2' 'Fz'}); % mpfc
chan=ismember(GA_TEPs_ACTIVE_ALL_PRE.label, {'P1' 'CP5' 'P3' 'CP3' 'PO3'}); % ag
gmfaAve.AVE = mean( dataMean(chan,:),1); 

% Plot gmfaAve.AVE to see which (trough and/or peak) looking for
figure, plot(time,gmfaAve.AVE)

for sidx = 1:size(site,1)
    % Find time points to define time series
    [val,tpW(1,1)] = min(abs(time-peakWin(1,1)));
    [val,tpW(1,2)] = min(abs(time-peakWin(1,2)));

    % Extract time series
    tseries = gmfaAve.(site{sidx,1});
    % Set -2 to 10 ms to 0;
    [~,tp1] = min(abs(-5-time));
    [~,tp2] = min(abs(17-time));
    tseries(1,tp1:tp2) = 0;

    latHold = [];
    num = 1;
    for b = tpW(1,1):tpW(1,2)       
        ztime = time(1,b);    
        % Find +/- windows
        tPlus = [];
        tMinus = [];
        if time(1,b) <= peakDefWin(1,2)
            for c = 1:peakDef(1,1)
                tPlus(c,1) = tseries(1,b) - tseries(1,b+c);
                tMinus(c,1) = tseries(1,b) - tseries(1,b-c);
            end
        else
            for c = 1:peakDef(1,2)
                tPlus(c,1) = tseries(1,b) - tseries(1,b+c);
                tMinus(c,1) = tseries(1,b) - tseries(1,b-c);
            end
        end

        %Find time points greater than 0 (peak)
        tPlusLog = tPlus > 0;
        tMinusLog = tMinus > 0;
        
%         %Find time points less than 0  (trough)
%         tPlusLog = tPlus < 0;
%         tMinusLog = tMinus < 0;
       
        testOut(1,b-(tpW(1,1)-1)) = size(tPlus,1) + size(tMinus,1);
        if  time(1,b) <= peakDefWin(1,2)
            if sum(tPlusLog) + sum(tMinusLog) == peakDef(1,1)*2
                latHold(num,1) = b;
                num = num+1;
            end
        else
            if sum(tPlusLog) + sum(tMinusLog) == peakDef(1,2)*2
                latHold(num,1) = b;
                num = num+1;
            end
        end
    end
    
    % Calculate latencies
    latenciesTemp = time(1,latHold);
end

% combine peaks
latencies.AVE = [latencies.AVE,latenciesTemp];  
latencies.AVE = sort(latencies.AVE);

% Calculate latency ranges and plot GMFAs
% Ranges are calcualted by taking half or the difference between subsequent
% peaks. First and last points are set to the time window of interest.
peakWin = peakWinOrig;
cond = condOrig;
gmfaAve = gmfaOrig;
for sidx = 1:size(cond,1)
    difference.(cond{sidx,1}) = diff(latencies.(cond{sidx,1}));
    
    for idx = 1:length(latencies.(cond{sidx,1}))
        if idx == 1
            latRange.(cond{sidx,1})(idx,1) = peakWin(1,1);
            latRange.(cond{sidx,1})(idx,2) = latencies.(cond{sidx,1})(1,idx)+ceil(round(difference.(cond{sidx,1})(1,idx))./2);
        elseif idx == length(latencies.(cond{sidx,1}))
            latRange.(cond{sidx,1})(idx,1) = latencies.(cond{sidx,1})(1,idx)-floor(round(difference.(cond{sidx,1})(1,idx-1))./2);
            latRange.(cond{sidx,1})(idx,2) = peakWin(1,2);
        else
            latRange.(cond{sidx,1})(idx,1) = latencies.(cond{sidx,1})(1,idx)-floor(round(difference.(cond{sidx,1})(1,idx-1))./2);
            latRange.(cond{sidx,1})(idx,2) = latencies.(cond{sidx,1})(1,idx)+ceil(round(difference.(cond{sidx,1})(1,idx))./2);
        end
    end
end

% Save latency ranges
save('peak_latency_ranges_ag.mat','latencies','latRange');
