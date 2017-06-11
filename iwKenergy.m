function ke = iwKenergy(u, v)
% ke = IWKENERGY(u, v)
%
%   inputs:
%       - u: x-component of the velocity.
%       - v: y-component of the velocity.
%
%   outputs:
%       - ke: kinetic energy density (J/m^2).
%
% IWKENERGY computes the horizontal kinetic energy density,
% rho0/2 * (u^2 + v^2).
% 
% Olavo Badaro Marques, 28/Nov/2016.


%% Define reference density value:

rho0 = 1025;


%% Compute kinetic energy density:

ke = (rho0/2) * (u.^2 + v.^2);

