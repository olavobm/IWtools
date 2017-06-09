function pe = iwPenergy(x, lrhoeta, n2backg)
% pe = IWPENERGY(x, lrhoeta, n2backg, z)
%
%   inputs:
%       - x: internal-wave driven density OR displacement perturbation.
%       - lrhoeta: TRUE (FALSE) if x is density (displacement)
%                  perturbation.
%       - n2backg: background buoyancy frequency squared field.
%
%   outputs:
%       - pe: potential energy.
%
% IWPENERGY computes the avilable potential energy density.
%
% Density perturbation (rho) and displacement (eta) can be easilly
% calculated from one another in a LINEAR framework as:
%               rho = (rho0/g) * N^2 * eta.
%
% Potential energy (in a LINEAR framework) is computed from rho or eta as:
%	PE = (0.5*g^2/rho0) * (rho/N)^2     or      PE = 0.5 *rho0*(N*eta)^2
%
% Displacement is usually computed from density perturbation, so this
% function computes APE from rho when lrhoeta is true. There may be
% some numerical details when computing rho or eta in different ways,
% which is the reason I included both options for computing APE in this
% function.
%
% Olavo Badaro Marques, 28/Nov/2016.


%% Define reference density value and gravitational acceleration:

rho0 = 1025;
g = 9.8;


%%

if ~isequal(size(n2backg), size(x))
    n2backg = repmat(n2backg, 1, size(x, 2), size(x, 3));
end


%% Compute potential energy density:

if lrhoeta
    pe = (0.5*g^2/rho0) * (x.^2)./n2backg;
else
    pe = (0.5*rho0) * (n2backg .* x.^2);
end

