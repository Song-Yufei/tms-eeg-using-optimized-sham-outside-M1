function [cleanSignal, iniExpInd, endExpInd, f0] = unDecaying(signal, timeToZP, ...
    baselineCorrTime, fsample, time, zeroSample, timeToFit, weightTime, tolTime, nChngTime, lowerLim, upperLim)
% Author Maryam Rostami. details please read the supplementary material: 1.1 Decay artifact removal
%___________________________________________
tCut1 = dsearchn(time', timeToZP(1));
tCut2 = dsearchn(time', timeToZP(2));
ind500 = dsearchn(time', nChngTime);
indEnd = length(time);

tolTimeSample = tolTime * fsample;
weightSample = weightTime * fsample;

signal(1:tCut1) = detrend(signal(1:tCut1),1); % detrend the baseline to avoid drops/jumps after tCut2
signal(tCut2:end) = detrend(signal(tCut2:end),1);
signal = signal - mean(signal(1+round(zeroSample+baselineCorrTime(1)*fsample):round(zeroSample+baselineCorrTime(2)*fsample)-1)); % baseline correction
signal(tCut1:tCut2) = signal(tCut2);

% >>>>>> exponential fit <<<<<<<<<
g = fittype('b*exp(c*x)'); % for fit function
cleanSignal = signal;
f0 = [];

maxTime = 1000e-3; % longest time point for decay trending search, No More used 
endExpPossibleInd = find(time>=maxTime,1); % previously cut1, the maximum ending point of the decay
signalSmoothed = smooth(signal(tCut2:end),20e-3*fsample);
slopeEst = diff(signalSmoothed);

if isempty(timeToFit)
    
    % figure; plot(time(tCut2:endExpPossibleInd),signal(tCut2:endExpPossibleInd),time(tCut2:endExpPossibleInd),signalSmoothed)    
    x11 = tCut1; y11 = signal(tCut1);
    x22 = tCut2+tolTimeSample; y22 = signal(tCut2+tolTimeSample);
    
    % use triangle method to find the start of the exp function (the initial exp index):
    % define point to line distance : https://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    pointsToLineDist = arrayfun(@(i) abs((x22-x11)*(y11-signal(i)) - (x11-i)*(y22-y11)) / sqrt((x22-x11)^2+(y22-y11)^2) ,x11:x22,'un',0);
    pointsToLineDist = [pointsToLineDist{:}];
    pointsAboveLine = arrayfun(@(i) (signal(i)-y11)/(i-x11) > (y22-y11)/(x22-x11), x11:x22, 'un', 0);
    pointsAboveLine = [pointsAboveLine{:}];

    endExpIndTemp = tCut2+50e-3*fsample; % this is only used for estimating the sign of exponential (rising/falling) and it's the first 50 msec. after tCut2
    slopeEstAvg = mean(diff(signalSmoothed(1:endExpIndTemp-tCut2+1)));
    if slopeEstAvg <= 0
        [~, maxDistIndTemp] = max(pointsToLineDist(pointsAboveLine));
        tempInds = find(pointsAboveLine);
    else
        pointsBelowLine = ~pointsAboveLine;
        [~, maxDistIndTemp] = max(pointsToLineDist(pointsBelowLine));
        tempInds = find(pointsBelowLine);
    end
    maxDistInd = tempInds(maxDistIndTemp);
    iniExpInd = maxDistInd-1+x11;
    if iniExpInd < tCut2
        iniExpInd = tCut2;
    elseif isempty(iniExpInd)
        iniExpInd = tCut2;
    end
    endExpInd = length(time); % 
else
    iniExpInd = timeToFit(1)*fsample + zeroSample;
    endExpInd = timeToFit(2)*fsample + zeroSample;    
end

%%%%%%%%%%%%%%%


if  ~isempty(iniExpInd) && ~isempty(endExpInd)  % do the fit if you found indices for the start and end of the exponential decays
    
    % fit function needs at least 3 samples (5 is conservative) and fit start should be earlier than tCut2+tolTimeMs
    % the "real" decay should be at least 5 msec.
    if (endExpInd-iniExpInd > 5e-3*fsample) && (iniExpInd <= tCut2+tolTimeSample) 

        cut11 = endExpInd - 1e-3*fsample;
        cut12 = endExpInd + 1e-3*fsample;
        zeroReach = [];
        yDiff_raw = abs(diff(signal([iniExpInd,endExpInd])));
        yDiff_smooth = abs(diff(signalSmoothed([iniExpInd-tCut2+1,endExpInd-tCut2+1])));
        if ~(yDiff_smooth<5) && ~(yDiff_raw<5) % if only the change in y axis is larger than 5 uV, then some exp decay might be seen, otherwise it's clean
                lenTemp = endExpInd-iniExpInd+1;
                if slopeEstAvg>0
                    f0 = fit(time(iniExpInd:endExpInd)', signal(iniExpInd:endExpInd)', g, ...
                        'TolFun',1e-4, 'TolX',1e-4,'StartPoint', [-200 -200],...%'StartPoint', [50 -200 -200],...
                        'Weights',[1e10 1e4*ones(1,weightSample-1) ones(1,lenTemp-weightSample-1) 1e4],...
                        'lower', lowerLim, 'upper', upperLim);

                else
                    f0 = fit(time(iniExpInd:endExpInd)', signal(iniExpInd:endExpInd)', g, ...
                        'TolFun',1e-4, 'TolX',1e-4,'StartPoint', [200 -200],...%'StartPoint', [50 -200 -200],...
                        'Weights',[1e10 1e4*ones(1,weightSample-1) ones(1,lenTemp-weightSample-1) 1e4],...
                        'lower', lowerLim, 'upper', upperLim);
                end
                zeroReach = endExpInd;
                if zeroReach > iniExpInd
                    yFitted = f0(time(iniExpInd:zeroReach))';
                    unDecayed = signal(iniExpInd:zeroReach) - yFitted;

                    lateSigblFixed = signal(zeroReach:end) - (signal(zeroReach)-unDecayed(end));
                    cleanSignal(zeroReach:end) = lateSigblFixed;

                    cleanSignal(iniExpInd:zeroReach) = unDecayed;    
                end              
        end
    end

    cleanSignal(tCut1:iniExpInd) = 0;
   
elseif isempty(iniExpInd)
    iniExpInd = NaN;

elseif isempty(endExpInd)
    endExpInd = NaN;
end

 cleanSignal(tCut1:iniExpInd) = 0;

% double check if the fit model was necessary and has not changed the data
% 1.considerably (mean of signal 500-end should not be changed)
% 2.fit model shouldn't be a straight line
% if abs((mean(cleanSignal(ind500:indEnd)) - mean(signal(ind500:indEnd)))) > 5 % set it to 5
%     cleanSignal = signal;
% end

if ~isempty(f0)
    xMid = time(round(mean([iniExpInd,endExpInd])));
    yMid = f0(xMid);
    slope1 = diff(f0(time([iniExpInd,endExpInd])))/diff(time([iniExpInd,endExpInd]));
    slope2 = diff([f0(time(iniExpInd))',yMid])/diff([time(iniExpInd),xMid]);
        
       if ( abs((mean(cleanSignal(ind500:indEnd)) - mean(signal(ind500:indEnd)))) > 5  && ...
            abs(mean(signal(iniExpInd:iniExpInd+20e-3*fsample)) - mean(f0(time(iniExpInd:iniExpInd+20e-3*fsample))))  > 2 ) || ...
            round(slope1)==round(slope2)  
        
            cleanSignal = signal;
            cleanSignal(tCut1:iniExpInd) = 0;

      end
end

