function conv_factor = factor_p2u(wvfreq, f, k)
% conv_factor = FACTOR_P2U(wvfreq, f, k)
%
%   inputs
%       - wvfreq:
%       - f:
%       - k:
%
%   outputs
%       - conv_factor:
%
% Conversion factor of the polarization relation
% between pressure and horizontal velocity (the component
% in the horizontal direction of wave propagation)
% (for normal modes only???).
%
% reference???
%
% Olavo Badaro Marques, 27/Jun/2018.


%%

rho0 = 1025;

conv_factor = (wvfreq .* k) ./ (rho0 .* (wvfreq.^2 - f.^2));