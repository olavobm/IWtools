function [mdsAmp, vmodes, fiterr] = fitVmodes(z, x, vmodes, zmds, xerror)
% [mdsAmp, vmodes, fiterr] = FITVMODES(z, x, vmodes, zmds, xerror)
% 
%   inputs
%       - z: depth grid vector associated with the data x.
%       - x: data whose columns are projected onto normal modes.
%       - vmodes: matrix whose columns are normal modes.
%       - zmds (optional): depth grid vector associated with the vertical
%                          modes, in case it is different than the data x.
%       - xerror (optional):
%
%   outputs
%       - mdsAmp: modal amplitudes array. One row for each mode (or
%                 column of vmodes) and one column for each of x.
%       - vmodes: return the matrix of vertical modes, which is useful if
%                 the vmodes input is not in the same z grid as the data.
%       - fiterr: struct variable containing quantities associated with
%                 the goodness of the fit. Its fields are:
%                   * xRes: residual of the fit
%                   * R2: variance explained by the fit.
%                   * merror: error on each model parameter.
%
% Function FITVMODES project the data x onto the normal modes specified
% by vmodes (i.e. each column of vmodes input represents a mode).
%
% If vmodes are available in a different depth grid than x, then
% specify zmds so the modes can be linearly interpolated onto the
% same z grid where is x given.
%
% Olavo Badaro Marques, 7/Nov/2016.


%%
if ~isvector(z)
    error(['Data depth array is NOT a vector. I could potentially ' ...
           'allow it to be a matrix, but it will require some work'])
end
    


%% Check whether a z-grid (zmds) was specified for the
% modes. If not, assumes it is the same as the data, but
% they must have the same number of rows:

if ~exist('zmds', 'var') || isempty(zmds)
    
    if size(vmodes, 1) == size(x, 1)
        zmds = z;
    else
        error(['Optional input zmds was not specified, ' ...
               'but number of rows of the modes is different than x'])
    end
end


%% Get the number of modes and columns of x:

nmds = size(vmodes, 2);
n = size(x, 2);


%% If necessary, interpolate vertical mode to the same depth
% level of the data to be projected onto the modes:

% check if modes have NaNs.....probably MUST not have!
% if modes do not go to the surface and bottom, interp1 below may give NaNs
% if it tries to extrapolate (in the case where min(zmds)>min(z) ||
% max(zmds)<max(z)

if ~isequal(z, zmds)
    
    vmodesinterp = NaN(length(z), nmds);

    for i = 1:nmds
        vmodesinterp(:, i) = interp1(zmds, vmodes(:, i), z);
    end
    
    vmodes = vmodesinterp;
    
end
   

%%

if exist('xerror', 'var')
    merror = NaN(nmds, n);
    lerror = true;
else
    lerror = false;
end

xRes = NaN(size(x, 1), n);
fitR2 = NaN(1, n);

mdsr2 = NaN(nmds, n);


%% Do the modal fit:

% Pre-allocate space for the modal amplitudes:
mdsAmp = NaN(nmds, n);

% Loop through the columns of x:
for i = 1:n
   
    % Check whether the ith column has enough
    % data points for the calculation:
    lgood = ~isnan(x(:, i));
    ngood = length(find(lgood));
    
    if ngood <= nmds
        warning(['Column i = ' num2str(i) ' does not have ' ...
                 'enough non-nan data points'])
        continue   % skip the fit and go to the next iteration of the loop
    end   
    
    % Project the ith column of x onto the normal modes:
    Gmodes = vmodes(lgood,  :);
    Gaux = (Gmodes' * Gmodes) \ (Gmodes');
    m = Gaux * x(lgood, i);
    
    % Assign m to ith column of mdsAmp:
    mdsAmp(:, i) = m;
    
    % Calculate the residue (difference between data and fit)
    xRes(lgood, i) = (Gmodes * m) - x(lgood, i);
    
% % %     %
% % %     fitR2(i) = 1 - (sum(xRes(lgood, i).^2) ./ ...
% % %                     sum((x(lgood, i) - mean(x(lgood, i))).^2));
    
	% Compute the correlation squared between each mode and the data.
    for i2 = 1:nmds
        mdsr2(i2, i) = corr(x(lgood, i), ...
                            mdsAmp(i2, i) .* Gmodes(:, i2)).^2;
    end
                
    % One could also use the misfit to calculate the error of the entire
    % fit (I'm not sure if that would work mode by mode)
    if lerror
        
        %
% %         xerror = sqrt(mean((xRes(lgood, i)).^2));
        
        %
        mCovmatrix = Gaux * ((xerror.^2) .* eye(ngood)) * Gaux';
        merror(:, i) = sqrt(diag(mCovmatrix));
    end

    
end

% Assign variables to error output structure
fiterr.misfit = xRes;
% % fiterr.R2 = fitR2;
fiterr.mdsr2 = mdsr2;
if lerror
    fiterr.merror = merror;
end


