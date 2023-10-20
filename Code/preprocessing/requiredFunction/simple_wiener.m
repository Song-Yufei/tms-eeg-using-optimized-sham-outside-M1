function [y_solved, sigmas] = simple_wiener(data,rf, sigmas_only)

% Estimating noise, undercorrelated over channels. See Mutanen, T. P., ...
% Metsomaa, J., Liljander, S., & Ilmoniemi, R. J. (2018). Automatic and ...
% robust noise suppression in EEG and MEG: The SOUND algorithm. ...
% Neuroimage, 166, 135-151.
%
% Input:
% data = (channels X times X trials) data matrix
% rf = regularization factor
% sigmas_only = 'true' if the only noise variance is needed
% 
% output:
% y_solved = filtered data of same dimensions
% sigmas = variances of noise for all channels
%

        [Nc, Nt, Nr]=size(data);
        data=reshape(data, Nc, []);
        C = cov(data');
        c=sum(diag(C));
        y_solved=zeros(Nc, Nt* Nr);
        if sigmas_only
            y_solved=[];
            sigmas=zeros(size(data,1),1);
             for i=1:Nc
                %disp(['Channel: ' num2str(i)])
                idiff = setdiff(1:Nc,i);
                sigmas(i) = C(i,i)-C(i,idiff)*((C(idiff,idiff)+c*rf/(Nc-1)*eye(Nc-1))\C(idiff,i));
            end
        else
        for i=1:Nc
            %disp(['Channel: ' num2str(i)])
            idiff = setdiff(1:Nc,i);
            y_solved(i,:) = C(i,idiff)*((C(idiff,idiff)+c*rf/(Nc-1)*eye(Nc-1))\data(idiff,:));
        end
        sigmas = (diag((data-y_solved)*(data-y_solved)'))/(size(data,2));
        y_solved=reshape(y_solved, [Nc, Nt, Nr]);
        end
        
        
