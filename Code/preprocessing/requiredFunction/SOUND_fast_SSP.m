function [corrected_data, sigmas,dn, correctionMatrix] = SOUND_fast_SSP(...
    data, LFM, iter,lambda0,U, inputData, sigmas, estimate_just_noise)

% This function performs simultaneous SSP_SIR and SOUND: SSP is performed
% within SOUND iterations. 

% INPUT:
%
% data: (channels x time points x trials) EEG
% LFM: lead-field matrix in the same reference as data
% iter: number of SOUND iterations
% lambda0: noise-to-signal ration, e.g., .01-.1 for averaged data, .5-10 for
% non-averaged
% U: out-projected signal-space dimensions in the columns (channels x
% dimensions)
% sigmas: known initial guess for the noise levels (variances) (channels x
% 1)
% estimate_just_noise: boolean-valued variable, if true, only sigmas are
% given as output
%
% OUPUT:
%
% corrected_data: same dimensions as input data
% sigmas: noise levels over channels (variances) (channels x 1)
% dn: relative change in noise levels in each iteration for convergence
% checking
% correctionMatrix: Spatial filter matrix which can be used to clean the
% original data by multiplying from the left
%
% Described in: Adapted beamforming: A robust and flexible
%approach for removing various types of artifacts from TMS–EEG data
%Johanna Metsomaa et. al. Submitted in 2023.
% .........................................................................
% 13 October 2023 : Johanna Metsomaa, Aalto university  
% .........................................................................

[n0, t0, r0]=size(data);
data=reshape(data, n0, []);

chanN = size(data,1);
if nargin < 7
    
sigmas = ones(chanN,1);
[~, sigmas] = simple_wiener(data,1e-4, true);
end
if nargin<4
    lambda0 = 1;
end

if nargin < 8
    estimate_just_noise = 0;
end

% Number of time points
T = size(data,2);



LL = LFM*LFM';

chanPerms = zeros(size(data,1),size(data,1)-1);
%Pall=zeros( n0-1, n0-1, n0);
Pin=eye(n0)-U*U';
for i = 1:size(data,1)
    is=setdiff(1:chanN,i);
    chanPerms(i,:) =is;
    
    [Us, ~, ~]=svds(Pin(is,:), n0-size(U,2)-1);
    P=Us';
    Pallin(:,:,i)=P;
    PLLP(:,:,i)=P*LL(is,is)*P';
    PL(:,:,i)=P*LFM(is,:);
    
    [Us, Ds, Vs]=svds(U(is,:), size(U,2));
    OperatorOut(:,:,i)=U(i,:)*Vs*diag(diag(Ds).^-1)*Us';
end

% Going through all the channels as many times as requested
dataCov=data*data'./size(data,2);

PLFM=Pin*LFM;
for k=1:iter
    sigmas_old = sigmas;
    
    %Evaluating each channel in a random order
    for i=1:n0
        chan = chanPerms(i,:);
            % Defining the whitening operator with the latest noise
            % estimates
            %W = diag(1./sigmas);
            
            % Computing the whitened version of the lead field
            %WL = (W(chan,chan))*(LFM(chan,:));
           % WL =  squeeze(Pall(:,:,i))*LFM(chan,:);
           % WLLW = W(chan,chan)*( LL(chan,chan)*(W(chan,chan))' );
            % WLLW =  *( LL(chan,chan)*(W(chan,chan))' );
             %lambdaCov=PLL0(chan,chan);
             lambda=lambda0*trace(squeeze(PLLP(:,:,i)))/(chanN-1);
            % wIn=(PLFM(i,:)*squeeze(PL(:,:,i))')*...
             %     ( (squeeze(PLLP(:,:,i)) + squeeze(Pallin(:,:,i))*lambda*diag(sigmas(chan))*squeeze(Pallin(:,:,i))')\squeeze(Pallin(:,:,i)));
              
             wIn=((squeeze(PLLP(:,:,i)) + ...
                 squeeze(Pallin(:,:,i))*lambda*diag(sigmas(chan))*squeeze(Pallin(:,:,i))')\...
                 (squeeze(PL(:,:,i))*PLFM(i,:)'))'*squeeze(Pallin(:,:,i)); 
              
             wOut=squeeze(OperatorOut(:,:,i));
             wM(i)=-1;
             wM(chan)=wIn+wOut;
            sigmas(i) =(wM*dataCov*wM');% sqrt((y_solved(i,:)-data(i,:))*(y_solved(i,:)-data(i,:))')/sqrt(T);
             
           
    % Following and storing the convergence of the algorithm
    dn(k) = max(abs(sigmas_old - sigmas)./sigmas_old);

    end
end

if estimate_just_noise
    x = [];
    corrected_data = [];

else
% Final data correction based on the final noise-covariance estimate.

          
            [Us, ~, ~]=svds(Pin, n0-size(U,2)-1);
            Pin=Us';
            correctionMatrix=((Pin*(LL + lambda0*trace(Pin*LL*Pin')/chanN*diag(sigmas))*Pin')\(Pin*LL))'*Pin;
            
            %corrected_data = (LFM*(Pin*LFM)')*((Pin*(LL + lambda0*trace(Pin*LL*Pin')/chanN*diag(sigmas))*Pin')\(Pin*data));
            
            if ~isempty(inputData)
                [n1, n2, n3]=size(inputData);
                corrected_data=reshape(correctionMatrix*reshape(inputData, n1, []), [n1, n2, n3]);
            else
                corrected_data=correctionMatrix*data;
            
                corrected_data=reshape(corrected_data, n0, t0, r0);
            end
            
            
end


end


