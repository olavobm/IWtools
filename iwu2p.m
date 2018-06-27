function p = iwu2p(u, wvfreq, f, k)
% p = IWU2P(u, wvfreq, f, k)
%
%   inputs
%       - u:
%       - wvfreq:
%       - f:
%       - k:
%
%   outputs
%       - p: pressure perturbation.
%
% The opposite of iwp2u.m
%
% See also: iwp2u.m
%
% Olavo Badaro Marques, 27/Jun/2018.

%%

p2u_factor = factor_p2u(wvfreq, f, k);

p = u ./ p2u_factor;
