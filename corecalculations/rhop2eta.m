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
n2 = repmat(n2, 1, size(rhop, 2));


%% Isopycnal displacement calculation:

eta = (rhop./n2) * (g/rho0);

