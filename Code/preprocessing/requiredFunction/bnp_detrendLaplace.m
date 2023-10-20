function [Xtr3d]=bnp_detrendLaplace(Xtr3,lambda)
%
% This function detrends the given EEG data epochs by Laplacian detrending
% See description in: Metsomaa, J, et al. "Causal decoding of individual 
% cortical excitability states." NeuroImage 245 (2021): 118652.
% 
% Input:
% Xtr3 = (channels X times X trials) EEG data
% lambda = regularization factor
%
% Output:
% Xtr3d = detrended data of the same size as input data
%
% .........................................................................
% 13 October 2023 : Johanna Metsomaa, Aalto university  
% .........................................................................


[Nc,Nt, Nr]=size(Xtr3);
I = speye(Nt);
D2 = spdiags(ones(Nt-2,1)*[1 -2 1],[0:2],Nt-2,Nt);

Xtr3d=zeros(Nc, Nt, Nr);
for i=1:Nr
trend = ((I+lambda^2*(D2'*D2))\squeeze(Xtr3(:,:,i))')';
Xtr3d(:,:,i)=Xtr3(:,:,i)-trend;

end

Xtr3d=Xtr3d-mean(Xtr3d,2);