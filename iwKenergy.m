function [ke, keint] = iwKenergy(u, v, z)
% [ke, keint] = IWKENERGY(u, v, z)
%
%   inputs:
%       - u: one horizontal velocity component.
%       - v: horizontal component perpendicular to u.
%       - z (optional): depth grid
%
%   outputs:
%       - ke: kinetic energy density.
%       - keint: depth-integrated kinetic energy density.
%
% IWKENERGY computes the horizontal kinetic energy density,
% rho0/2 * (u^2 + v^2). If optional input z is provided, IWKENERGY
% also computes the depth-integrated kinetic energy density.
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


%% Compute depth-integrated quantity:

% Only compute if input z was specified:
if exist('z', 'var')
    
    % If u or v has NaNs, prints a warning message. Create
    % logical value to tell whether there are NaNs:
    if ~isempty(find(isnan(u), 1)) || ~isempty(find(isnan(v), 1))
        warning(['Velocity components have NaN. Depth-integrated ' ...
                 'quantity may be very different than true value.'])
             
        lnan = true;
    else
        lnan = false;
    end
    
    % If there are no NaNs, integrate right away:
    if ~lnan
        keint = trapz(z, ke, 1);
        
	% If there are NaNs, integrate each column separately:
    else
        ncols = size(ke, 2);

        keint = NaN(1, ncols);
        for i = 1:ncols
            lok = ~isnan(ke(:, i));
            keint(i) = trapz(z(lok), ke(lok, i), 1);
        end
    end
    
end