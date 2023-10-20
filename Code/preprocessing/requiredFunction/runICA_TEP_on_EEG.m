function [EEG_out] = runICA_TEP_on_EEG(EEG_in, ts_out)

% This function wraps up the ICA code provided by Johanna Metsomaa. 

mixedsig = EEG_in.data;

ts = EEG_in.times;
ts_in = [min(ts),max(ts)];
if isfield(EEG_in,'tmscut') 
    ts_out = EEG_in.tmscut.cutTimesTMS;
    [~,t0ind] = min(abs(EEG_in.times-ts_out(1)));
    [~,t1ind] = min(abs(EEG_in.times-ts_out(2)));
elseif isempty(ts_out)
    ts_out = [];
end

chanlocs = EEG_in.chanlocs;
maxNumIter = 1000;
rmch = zeros(1,length(chanlocs));
tempMean = [];

[A, indsRem, icasig, icasphere, W] = runICA_TEP(mixedsig,ts, ts_in, ...
    ts_out, chanlocs, maxNumIter, rmch, tempMean);

%icasig = cat(2, cat( 2, icasig(:,1:(t0ind-1),:), zeros(size(icasig,1), t1ind - t0ind +1, size(icasig,3) ) ),icasig(:,(t1ind+1):end,:) );
if isempty(ts_out)

else
    icasig = [icasig(:,1:(t0ind-1),:), zeros(size(icasig,1), t1ind - t0ind +1, size(icasig,3) ) ,icasig(:,(t0ind):end,:) ];
end
disp('after correctin')
size(icasig)
% % Sorting the components based on the variance:
% 
% for x = 1:size(icasig,1)
%     vars(x) = var(mean(icasig(x,:,:),3));
% end
% varsPerc = vars/sum(vars)*100;
% 
% [~, sortedICs] = sort(varsPerc,'descend');


EEG_out = EEG_in;
%EEG_out.data = x_corr; % not needed, because nothing is removed
EEG_out.ICA_mixing_matrix = A; % mixing matrix
EEG_out.removed_ICs = indsRem; 
EEG_out.icaact=icasig; % EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
EEG_out.icawinv=A;
EEG_out.icasphere=icasphere;
EEG_out.icaweights=W;
EEG_out.icachansind=1:length(chanlocs);

% Ranks and sorts components
% EEG_out = tesa_sortcomps(EEG_out);

end

