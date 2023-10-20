function [Xdt]=detrendAllData(Xorig, ts, n, excludeInterval)
%
% performs 2-rounds of robust detrending of TEPs (or other epoched data)
% over all channels
% 1st round is to find the offset level robustly
% 2nd round is to find the polynomial of given order
%
% input:
% Xorig: original data channels x times x trials
% ts: time axis 
% n: order of the fitted polynomial, e.g. 3
% excludeInterval: time interval including the evoked response, having ...
% non-zero average signal, which should be excluded from the fit
% [start_time, end_time], e.g. [-4, 600]
% If your data do not have evoked activity set end_time < start_time, e.g.,
% [-1 0]
% 
% output:
% Xdt: detrended data channels channels x times x trials
%
% Description in Hernandez-Pavon, Julio C., et al. "Removing artifacts ...
% from TMS-evoked EEG: A methods review and a unifying theoretical ...
% framework." Journal of Neuroscience Methods 376 (2022): 109591.
% .........................................................................
% 29 March 2021 : Johanna Metsomaa, BNP, University of Tuebingen  
% .........................................................................
[C, T, R]=size(Xorig);

Xdt=zeros(C,T,R);
for i=1:C
    
        [yestim, ~]=removePolyTrendlineTEP_robust(permute(Xorig(i,:,1:end),[2 3 1]), excludeInterval(1), excludeInterval(2), ts,n, false);
        Xdt(i,:,:)=yestim;
end
