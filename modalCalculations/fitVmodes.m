function [mdsAmp, vmodes, xRes] = fitVmodes(z, x, vmodes, zmds)
% [mdsAmp, vmodes, xRes] = FITVMODES(z, x, vmodes, zmds)
% 
%  inputs:
%    - z: depth grid vector associated with the data x.
%    - x: data whose columns are projected onto normal modes.
%    - vmodes: matrix whose columns are normal modes.
%    - zmds (optional): depth grid vector associated with the vertical
%                       modes, in case it is different than the data x.
%
%  outputs:
%    - mdsAmp: modal amplitudes array. One row for each mode (or column of
%              vmodes) and one column for each of x.
%    - vmodes: return the matrix of vertical modes, which is useful if the
%              vmodes input is not in the same z grid as the data.
%    - xRes:
%
% Function FITVMODES project the data x onto the normal modes specified
% by vmodes(i.e. each column of vmodes input represents a mode).
%
% If vmodes are available in a different depth grid than x (i.e. if
% observations of N2 are at different depths than the data x), please
% specify zmds so the modes can be linearly interpolated to the z grid
% associated with x.
%
% Olavo Badaro Marques, 7/Nov/2016.


%% Check whether a z-grid (zmds) was specified for the
% modes. If not, assumes it is the same as the data, but
% they must have the same number of rows:

if ~exist('zmds', 'var')
    
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
    m = ( vmodes(lgood,  :)' * vmodes(lgood,  :)) \ ...
        ( vmodes(lgood,  :)' * x(lgood, i));
    
    % Assign m to ith column of mdsAmp:
    mdsAmp(:, i) = m;
    
    
end




