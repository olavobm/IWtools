function eta = rhop2eta(rhop, N, zRhoN, rho0)
% eta = RHOP2ETA(rhop, N, zRhoN, rho0)
%
%   inputs:
%       - rhop: density perturbation (kg per m^3).
%       - N: time-averaged buoyancy frequency (radians per s).
%       - zRhoN (optional): 1x2 cell array.
%       - rho0 (optional): reference density.
%
%   outputs:
%       - eta: isopycnal displacement in meters.
%
% Based on a linear approximation, compute isopycnal
% displacement from a density anomaly (rhop) and a
% background buoyancy frequency profile (N).
%
% Olavo Badaro Marques, 28/Apr/2017.


%% Parameters:

g = 9.8;

if ~exist('rho0', 'var')
    rho0 = 1025;
end


%%

if ~exist('zRhoN', 'var') || isempty(zRhoN)
    
    NatRho = N;
    
else
    
    NatRho = interp1(zRhoN{2}, N, zRhoN{1});
    
end

n2 = NatRho.^2;


%% If the sizes of n2 and rhop are not the same, then
% match the size of n2 to be the same as rhop

n2size = size(n2);
rpsize = size(rhop);

% If the variables do NOT have the
% same size, than match n2 to rhop
if ~isequal(n2size, rpsize)
    
    if isvector(n2)
        
        n2 = n2(:);
        
        repdims = [1, rpsize(2:end)];
        
        n2 = repmat(n2, repdims);
        
    else
        
        % in this case, it may be that N2 is a
        % matrix and rhop is a 3-D array.
        
    end
        
end


%% Isopycnal displacement calculation:

eta = (rhop./n2) * (g/rho0);

