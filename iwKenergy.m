function ke = iwKenergy(u, v, rho0)
% ke = IWKENERGY(u, v, rho0)
%
%   inputs
%       - u: x-component of the velocity.
%       - v:   y- "       "  "     ".
%       - rho0 (optional): density reference.
%
%   outputs
%       - ke: kinetic energy density (J/m^2).
%
% IWKENERGY computes the horizontal kinetic energy density,
% rho0/2 * (u^2 + v^2).
% 
% Olavo Badaro Marques, 28/Nov/2016.


%% Define reference density value

if ~exist('rho', 'var')
    rho0 = 1025;
end


%% Compute kinetic energy density

ke = (rho0/2) * (u.^2 + v.^2);

