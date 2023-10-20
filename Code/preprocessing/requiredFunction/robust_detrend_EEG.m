function [EEG_out] = robust_detrend_EEG(EEG_in, polynomial_order)

% This function wraps up the robust detrend code provided by Johanna Metsomaa. 
% Tuomas Mutanen

EEG_out = EEG_in;
data_in = EEG_in.data;

n = polynomial_order;

size(data_in)

excludeInterval = [-20 600];


[data_out] = detrendAllData(data_in, EEG_in.times, n, excludeInterval);


EEG_out.data = data_out;

end

