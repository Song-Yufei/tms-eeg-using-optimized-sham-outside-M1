function [Sestim, SestimPow, times]=estimateCorrelatedSourceAmplitudesFast(LnrN, Mdw, Xst, lambda, indsp, timeAxis, timeInt)

% This function estimates source amplitudes by MNE of correlated source
%
% INPUT:
% LnrN - lead-field matrix (fixed orientations only in this function)
% channels X sources
% Mdw - dewhitening matrix, meaning Mdw*Mwd = source covariance matrix,
% (sources X sources)
% Xst - data of interest, (channels X times) 
% lambda - regularization factor
% indsp -indices of the included data (based on prior knowledge)
% timeAxis - the time for each data sample, 1 X time points (can also be 1:number_of_time_points 
% if this irrelevant)
% timeInt - [first_time_label, last_time_label] defining the period of
% interest (can also be, eg.,[1, number_of_time_points])
%
% OUTPUT:
% Sestim - source amplitudes
% SestimPow - powers of sources
% times - time axis for sources
% .........................................................................
% 29 March 2021 : Johanna Metsomaa, BNP, University of Tuebingen  
% .........................................................................

nOrigSources=size(LnrN,2);

if ~isempty(Mdw)
    
    LnrN=reshape(LnrN(:,indsp)*Mdw, [size(LnrN,1), length(indsp)]);   
else
    Mdw=eye(length(indsp));
end


[~, it1]=min(abs(timeAxis-timeInt(1)));
[~, it2]=min(abs(timeAxis-timeInt(2)));
times=timeAxis(it1:it2);

    Sestim=zeros(nOrigSources, it2-it1+1); %varataan muistia lähde-estimaateille
    SestimPow=zeros(nOrigSources, it2-it1+1); %varataan muistia lähde-estimaateille


%for i=1:(it2-it1+1)
    
    %disp(['Sample index: ' num2str(i)])
    %[B,~]=lassoglm(LnrN,Xst(:,it1+i-1),'normal','alpha', alpha, 'standardize', false, 'lambda', trace(Xst*Xst')/numel(Xst)*lambda); % [5]*1e1 lambda
    Sestim(indsp,:)=...
        ((Mdw*LnrN')/(LnrN*LnrN'+eye(size(LnrN,1))*trace(Xst*Xst')...
        /numel(Xst)*lambda))*Xst(:, it1:it2);

    SestimPow(indsp,:)=Sestim(indsp,:).^2;
        
 %end       

