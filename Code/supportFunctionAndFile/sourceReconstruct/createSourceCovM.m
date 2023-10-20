function [Mdw,indsp]=createSourceCovM(headmodel, distanceTH, attenuation_length)
% This function estimates correlation matrix (in spatial dimension) for
% source space based on L2 distances of the sources. It also eliminates
% deep sources, and computes spatial whitening matrix for the source space
% to generate a set of uncorrelated surrogate sources
%
% INPUT:
% headmodel - headmodel struct
% distanceTH - threshold, in cm, for choosing sources based on distance from the
% origin. Sources are included if they excee the th (try eg 6),  
% [the recommended distance for dipole sources from the inner skull is half
% of triangle side length (Matti hbf toolbox)]
% attenuation_length - (e.g. 15) the length when the cross-covariance
% decreases as defined by std in the Normal distribution
%
% output: 
% Mdw - whitening matrix used in the source estimation
% indsp - the kept source indices
% .........................................................................
% 29 March 2021 : Johanna Metsomaa, BNP, University of Tuebingen  
% .........................................................................

p=permute(headmodel.smesh.p, [1 3 2]); % source locations (x,y,z)
p=repmat(p,[1 size(p,1), 1]);
dists=(sum((p-permute(p,[2 1 3])).^2,3)); %distance between all source locations
SigmaB=exp(-dists./attenuation_length^2); %source covariance matrix with covariances decreasing by the Gaussian function

pd=sqrt(sum(headmodel.smesh.p.^2,2));
indsp=find(pd> distanceTH); % choose source with at least a specified distance from the orginin
[u, d, ~]=svd(SigmaB(indsp, indsp));
d=diag(d);
ncomp=length(indsp);
%Mw=diag(1./sqrt(d(1:ncomp)))*u(:,1:ncomp)'; 
Mdw=u(:,1:ncomp)*diag(sqrt(d(1:ncomp))); 
