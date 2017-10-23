function [xMode, xErr] = modalRecons(x, mdsAmp, vmodes)
%
%
%   inputs:
%       - x:
%       - mdsAmp:
%       - vmodes:
%
%   outputs:
%       - xModes:
%       - xRes:
%
% Olavo Badaro Marques, 06/June/2017.

% ------------------------------------------------------------------------
% variance explained by the modes could be either in time or depth....
% ------------------------------------------------------------------------


%%

% TO DO: should check if inputs sizes are consistent...


%%

% TO DO: extra input to subset the columns of vmodes...


%% If there are NaNs, the calculation might be optimized in the
% case where the NaNs are found in the same rows for all columns:

if any(isnan(x(:)))
    
    colsets = sepColsNanPos(bla);
     
    if length(colsets)==1
        ntypecols = 1;
%         indcols = 1:size(x, 2);
    else
        ntypecols = size(x, 2);
    end
    
end


%%

xMode = NaN(size(x));


%%

% Constant factor which determines the confidence interval (CI)
% associated with the error. 1.96 gives a 95% CI:
facCI = 1.96;

%
for i = 1:ntypecols
    
    %%
    
    %
    if ntypecols==1
        indcols = 1:1:size(x, 2);
    else
        indcols = i;
    end
    
    %
    lok = ~isnan(vmodes(indcols(1), :));
    
    %
    vmodesOK = vmodes(lok, :);
    
    xMode(:, indcols) =  * mdsAmp;
    
    
    %%
    
    %
    xErr.res = xMode - x;
    xErr.MSE = mean(xErr.res.^2);

    % Compute r2 (squared correlation coefficient):
%     xErr.r2 = var(xMode) / var(xgd);   % which dimension?????

    %
    errorVar = (facCI^2) * xErr.MSE;
    
    auxGmodes = (vmodesOK' * vmodesOK) \ vmodesOK';
    
    %
    xErr.mErr = errorVar .* diag(auxGmodes * auxGmodes');
    xErr.mErr = sqrt(err.mErr);
    
    
end



%% -------------------------------------------------------

end


%% -------------------------------------------------------
% --------------------------------------------------------
% --------------------------------------------------------

% function [] = reconsModes(G, m)
% 
% 
% end

