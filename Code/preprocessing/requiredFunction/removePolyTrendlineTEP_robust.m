function [yestim, ybaseline]=removePolyTrendlineTEP_robust(y, t1, t2, ts,n0, visual)
% performs 2-rounds of robust detrending of TEPs (or other epoched data)
% for one channel
% 1st round is to find the offset level robustly
% 2nd round is to find the polynomial of given order
%
% output:
% yestim: estimtated detrended signal (times x trials)
% ybaseline: the trend (times x trials)
% 
% input:
% y: the data from one channel (times x trials)
% t1: the start time of the artifact/responses
% t2: end time of the response (where the original baseline level has been
% reached). If your data do not have evoked activity set t2<t1, e.g., t1=0,
% t2=-1;
% ts: time axis (1 x times)
% n: order of the polynomial, e.g., 3. 
% visual: logical true/false if you want to visualize or not, respectively
%
% .........................................................................
% 29 March 2021 : Johanna Metsomaa, BNP, University of Tuebingen  
% .........................................................................


[~, t1]=min(abs(ts-t1));
[~, t2]=min(abs(ts-t2));


blC=mean(y(1:t1,:),1); %baseline level
y=y-repmat(blC,[size(y,1),1]); %remove offset

indsR=t1:(t2-1);%removed interval
ymask0=true(size(y));

ymask0(indsR,:)=false; %remove interval from mask0
ymask=ymask0;
%set starting variables for iterations
yestim=zeros(size(y));
ybaseline=yestim;
ybaselineOld=ybaseline;
continF=true;
in=1;
nall=[0 n0];
multipAll=[3 3];

witer=1;
inds=1:length(ts);
count=0;
cont3=true;
ytemp=y;

while continF
    n=nall(in);
    multipTH=multipAll(in);
    count=count+1;
    
for k=1:size(y,2)
    
    indsP=inds(ymask(:,k));
    
%     p2=polyfit(indsP,ytemp(ymask(:,k),k)', n);%% ???? transpose
    p2=polyfit(indsP',ytemp(ymask(:,k),k), n);%% Johanna's standard
%     scripts. the outputs are the same to the last code
    
    ybaseline( :, k) = polyval(p2,1:length(ts));
   
end


yestim=y-ybaseline;
yestimD=yestim(:);%(ymask);
yestimD=sort(yestimD);

nD=round(length(yestimD)*.01);
yestimD=yestimD((nD+1):(end-nD));
th=sqrt(mean(yestimD.^2))*multipTH;
if visual
subplot(2,1,1)
hold off
hist(yestimD)
xlabel('Difference to baseline')
ylabel('Count')
title(['Threshold: ' num2str(th) ' at iteration:  ' num2str(witer)])
subplot(2,1,2)
hold off
plot(ts,y(:,64)), hold on, plot(ts, ybaseline(:,64), 'k','linewidth', 1.5)
xlabel('Time (ms)')
ylabel('Amplitude (\muV)')
%size(ymask0), size(y), size(abs(yestim)>th)

end
ymask=true(size(y));
ymask(abs(yestim)>th)=false;


Diff=max(abs(ybaseline(:)-ybaselineOld(:)));

continF=((Diff>1 || cont3) || n<n0) && count<100;

% if Diff<1 && ~cont3
%     in=in+1;
%     cont3=true;
% end

if (Diff<1 || ~cont3)
    
if  n==n0
    
    ymaskN3=true(size(y));
    ymaskN1=true(size(y));
    
    ymaskN3(abs(yestim)>std(yestimD)*multipTH)=false;
    ymaskN1(abs(yestim)>std(yestimD)*1)=false;
for ir=1:size(y,2)
    v1=~ymaskN1(:,ir);
    v3=~ymaskN3(:,ir);
    nl = length(v1);
    e = ones(nl,1);
    Mn = spdiags([e e e],-1:1,nl,nl);
    Mn(~v1,:)=0;
    Mn(:,~v1)=0;
    Mn = Mn+spdiags([e],0,nl,nl);
 
   ymask(:,ir)=~logical((Mn)\v3);
   
end
    ymask(~ymask0)=false;
    cont3=false;    
else
    in=in+1;
    cont3=false;
end
end
if visual
    bar(ts(~ymask(:,64)), y(~ymask(:,64),64), 'edgecolor','r')
set(gca, 'ylim', [-70 70], 'xlim', [ts(1) ts(end)])
ginput(1);
end


ybaselineOld=ybaseline;
witer=witer+1;
%y=y-repmat(mean(ybaseline(:,1:t1,:),2),[1, size(y,2),1]);

end

yestim =y-ybaseline;
ybaseline=ybaseline+repmat(blC,[size(y,1),1]);





