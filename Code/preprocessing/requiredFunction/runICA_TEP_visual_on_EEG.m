function [EEG_out] = runICA_TEP_visual_on_EEG(EEG_in, ts_out, keepAll)

% This function wraps up the ICA code provided by Johanna Metsomaa. 
EEG_in.times(end)
mixedsig = EEG_in.data;

ts = EEG_in.times;
ts_in = [min(ts),max(ts)];
if isfield(EEG_in,'tmscut') 
    ts_out = EEG_in.tmscut.cutTimesTMS;
elseif isempty(ts_out)
    ts_out = [];
end

chanlocs = EEG_in.chanlocs;
maxNumIter = 10;
rmch = zeros(1,length(chanlocs));
tempMean = [];

[x_corr, ~ , A, indsRem, icasig, icasphere, W] = runICA_TEP_visual(mixedsig,ts, ts_in, ...
    ts_out, chanlocs, maxNumIter, rmch, tempMean, keepAll);

EEG_out = EEG_in;
EEG_out.times(end)
EEG_out.data = x_corr;
EEG_out.ICA_mixing_matrix = A;
EEG_out.removed_ICs = indsRem; 
EEG_out.icaact=icasig; % EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
EEG_out.icawinv=A;
EEG_out.icasphere=icasphere;
EEG_out.times(end)
EEG_out.icaweights=W;
EEG_out.icachansind=1:length(chanlocs);
end

