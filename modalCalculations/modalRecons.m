function [xMode, xErr] = modalRecons(x, mdsAmp, vmodes)
%
%
%   inputs:
%
%
%
%
%   outputs:
%       - xModes:
%       - xRes:
%
% Olavo Badaro Marques, 06/June/2017.


%% If there are NaNs, the calculation might be optimized in the
% case where the NaNs are found in the same rows for all columns:

if any(isnan(x(:)))
    
    colsets = sepColsNanPos(bla);
     
    if length(colsets)==1
        ntypecols = 1;
    else
        ntypecols = size(x, 2);
    end
    
end


%%





end


%% -------------------------------------------------------
% --------------------------------------------------------
% --------------------------------------------------------

function [] = reconsModes(G, m)


end



%%
% function [mdsAmp, vmodes, xRes] = fitVmodes(z, x, vmodes, zmds)
% -------------------------------------
    % Compute the misfit and the mean-squared error (MSE):
    fit4err = (Gmodes * m);
    err.res = fit4err - xgd;
    err.MSE = mean(err.res.^2);
    
    % Compute r2 (squared correlation coefficient):
    err.r2 = var(fit4err) / var(xgd);
    
    % Constant factor which determines the confidence interval (CI)
    % associated with the error. 1.96 gives a 95% CI. In fact, this
    % constant is valid if the variable is normally distributed:
    facCI = 1.96;
    
    % Use the MSE as the magnitude of data error variance. Define
    % the error variance and the matrix that combines with the data
    % error to give error estimates for the model parameters. Note
    % that, in general, an error covariance matrix should be used instead
    % of a scalar. But it is assumed the error is constant and
    % uncorrelated at different locations, which greatly decreases the
    % computational time:
    errorVar = (facCI^2) * err.MSE;
    aux_G4err = (G4err'*G4err) \ G4err';
    
    % Compute the error variance of the model parameters and take
    % the square root to compute the error associated with the
    % confidence interval set by facCI:
    err.mErr = errorVar .* diag(aux_G4err * aux_G4err');
    err.mErr = sqrt(err.mErr);