% This function estimate the subspace of artifacts (dependent on
% dataset provided)


function [Pexp, Uexp]=ProjectionMatrixComputation(Xnew, times,n, intTfit, chanlocs, Xave)
% Input:
% Xnew = (channels X times X trials) data by which PCA is used to estimate
% artifact subspace
% times = (1 x times) time axis
% n = use empty matrix [] as input here
% intTfit = (1 X 2) vector with the entries giving the start and end time
% point for plotting
% chanlocs = chanlocs struct from the EEGlab data structure for plotting
% Xave = (channels X times) averaged data over the epochs for visualizing
% the averaged component waveforms
% 
% Output:
% Uexp = artifact subspace
% Pexp = projection matrix to completely wipe out the data in this artifact
% subpace
% .........................................................................
% 13 October 2023 : Johanna Metsomaa, Aalto university  
% .........................................................................



Nc=size(Xnew,1);
[~,it1]=min(abs(times-intTfit(1)));
[~,it2]=min(abs(times-intTfit(2)));

if isempty(n)
    Xcov=reshape(Xnew, Nc,[]);
    [u, d, ~]=svd(Xcov, 'econ');
    v=u'*mean(Xave,3);
    
    b=1;
   figure('units','normalized','outerposition',[0 0 1 1])
        subplot(4,2,1)
    bar(diag(d))
    subplot(4,2,2)
    plot(times(it1:it2), Xave(:, it1:it2))
     [x,~]=ginput(1);
    while b==1
        
   
    n=round(x);
    for i2=1:n
    
        subplot(4,n,i2+n)
        topoplot(u(:,i2), chanlocs) 
        
        subplot(4,n,i2+n*2)
        plot(times(it1:it2), v(i2,it1:it2))
        
        Uexp=u(:,1:i2);
        Pexp=eye(Nc)-Uexp*Uexp';
        subplot(4,n,i2+n*3)
        plot(times(it1:it2), Pexp*Xave(:, it1:it2))
        
    end
    [x,~,b]=ginput(1);
    close gcf
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(4,2,1)
    bar(diag(d))
    subplot(4,2,2)
    plot(times(it1:it2), Xave(:, it1:it2))
    
    end
    n=round(x);
    Uexp=u(:,1:n);
    Pexp=eye(Nc)-Uexp*Uexp';
end
close gcf