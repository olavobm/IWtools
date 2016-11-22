function [mdsAmp, xRes] = fitVmodes(z, x, vmodes, zmds)
% = FITVMODES(z, x, vmodes, zmds)
% 
%  inputs:
%    - z: 
%    - x:
%    - vmodes:
%    - zmds (optional):
%
%  outputs:
%    - mdsAmp:
%    - xRes:
%
% I MUST RETURN THE VMODES IN THE GRID Z (I ALREADY COMPUTE IN THE CODE AND
% IT IS REQUIRED FOR FOR FURTHER CALCULATIONS WITH mdsAmp).
%
% Olavo Badaro Marques, 7/Nov/2016.


%%

if ~exist('zmds', 'var')
    zmds = z;
end


%% 

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
    

%% Convert the matrix vmodes into a 1xnmds cell array,
% which is the appropriate format 

imf.modelinput = mat2cell(vmodes, length(z), ones(1, nmds));


%%

% indnan = find(isnan(x));
% [inan, jnan] = ind2sub(size(x), indnan);

mdsAmp = NaN(nmds, n);


for i = 1:n
   
%     clear imf   % just to make sure
    
    lgood = ~isnan(x(:, i));
    ngood = length(find(lgood));
    
    if ngood <= nmds
        warning(['Column i = ' num2str(i) ' does not have ' ...
                 'enough non-nan data points'])
        continue   % skip the fit and go to the next iteration of the loop
    end
    
    % My least squares function already deals with NaNs:
%     imf.modelinput = mat2cell(vmodes(lgood, :), ngood, ones(1, nmds));
    [~, m] = myleastsqrs(z, x(:, i), imf);
    
    % Assign m to ith column of mdsAmp:
    mdsAmp(:, i) = m;
    
    
end




