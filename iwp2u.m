function u = iwp2u(p, wvfreq, f, k)
% u = IWP2U(p, wvfreq, f, k)
%
%   inputs
%       - p: pressure perturbation.
%       - wvfreq:
%       - f:
%       - k:
%
%   outputs
%       - u: 
%
% k is the horizontal wavenumber.
%
% For constant stratification, k of the normal modes is
% k = s*m, where m is the vertical wavenumber of the a
% normal mode and s is the internal wave characteristic
% slope.
%
%
% Olavo Badaro Marques, 27/Jun/2018.


%%

p2u_factor = factor_p2u(wvfreq, f, k);

u = p .* p2u_factor;