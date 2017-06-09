function ke = iwKenergy(u, v)
% ke = IWKENERGY(u, v)
%
%   inputs:
%       - u: one horizontal velocity component.
%       - v: horizontal component perpendicular to u.
%
%   outputs:
%       - ke: kinetic energy density.
%
% IWKENERGY computes the horizontal kinetic energy density,
% rho0/2 * (u^2 + v^2).
%
% If v input is not specified (such as in a synthetic case where
% the wave is aligned the axis), create a v array with zeros only.
% 
% Olavo Badaro Marques, 28/Nov/2016.


%% Define reference density value:

rho0 = 1025;


%% If v input was not specified or given as empty array, create
% a v array with zeros (that do not change the kinetic energy)

if ~exist('v', 'var') || isempty(v) 
    v = zeros(size(u));
end


%% Compute kinetic energy density:

ke = (rho0/2) * (u.^2 + v.^2);

