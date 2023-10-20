function [data_sound_filt, W]=FilterOutGivenPatternsSound_BF(A, data_sound, Msound)

% This function removes the artifact patterns, given as topographies in A,
% by beamforming-based spatial filterin. It is to be used for SOUND-cleaned
% data.
%
% Input:
% A = (channels X patterns) mixing matrix, where columns are the artifact
% topographies
% data_sound = (channels X times X trials) EEG data, which have been cleaned by
% SOUND 
% Msound = (channels X channels) cleaning matrix from the SOUND algorithm
%
% Output:
% W = (patterns X channels) beamforming filter matrix, which retrieves the
% artifact waveforms
% data_sound_filt = cleaned data after removing artifacts defined by A
%
% .........................................................................
% 13 October 2023 : Johanna Metsomaa, Aalto university  
% .........................................................................


[C, T, R]=size(data_sound);

cov=reshape(data_sound, C, []);
cov=cov*cov'/(T*R);
[u, d, ~]=svd(cov);
d=sqrt(diag(d));

figure,

bar(cumsum(d)/sum(d))
set(gca, 'xlim', [0 50], 'ytick', 0:.05:1) % [0 30] is original setting
hold on
plot([0 50], [.9 .9 ], 'r--', 'linewidth', 2) % threshold can be set manully, [0.85 0.85]
% the threshold determines how many components (PCA)required to project
% back to the channel level.
% if it is too small, then real frontal activity can be suppressed.
% [0 30] is original setting
[xn, ~, ~]=ginput(1);
xn=round(xn);
P=u(:,1:xn)';

icov=diag(1./d(1:xn));
pA=P*Msound*A;

W=inv(pA'*icov*pA)*pA'*icov*P;

data_sound_filt=data_sound-reshape((P'*pA)*(W*reshape(data_sound, C, [])),...
    C, T, R);
