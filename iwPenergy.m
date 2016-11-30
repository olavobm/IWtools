function [pe, peint] = iwPenergy(x, lrhoeta, n2backg, z)
% [pe, peint] = IWPENERGY(rho, n2backg)
%
%   inputs:
%       - x: internal-wave driven density OR displacement perturbation.
%       - lrhoeta: TRUE (FALSE) if x is density (displacement)
%                  perturbation.
%       - n2backg: background buoyancy frequency squared field.
%       - z (optional): depth grid
%
%   outputs:
%       - pe: kinetic energy.
%       - peint: depth-integrated kinetic energy.
%
% IWPENERGY computes the avilable potential energy density. If optional
% input z is provided, IWPENERGY also computes the depth-integrated
% available potential energy density.
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

%% Compute potential energy density:

if lrhoeta
    pe = (0.5*g^2/rho0) * (x.^2)./n2backg;
else
    pe = (0.5*rho0) * (n2backg .* x.^2);
end


%% Compute depth-integrated quantity:

% Only compute if input z was specified:
if exist('z', 'var')
    
    % If x has NaNs, prints a warning message. Create
    % logical value to tell whether there are NaNs:
    if ~isempty(find(isnan(x), 1))
        warning(['Input x has NaN. Depth-integrated ' ...
                 'quantity may be very different than true value.'])
             
        lnan = true;
    else
        lnan = false;
    end
    
    % If there are no NaNs, integrate right away:
    if ~lnan
        peint = trapz(z, pe, 1);
        
	% If there are NaNs, integrate each column separately:
    else
        ncols = size(pe, 2);

        peint = NaN(1, ncols);
        for i = 1:ncols
            lok = ~isnan(pe(:, i));
            peint(i) = trapz(z(lok), pe(lok, i), 1);
        end
    end
    
end