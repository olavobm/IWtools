function [KE, KEphasor] = phasors2KE(u, v)
% [KE, KEphasor] = PHASORS2KE(u, v)
%
%   inputs
%       -
%       -
%
%   inputs
%       -
%       -
%
%
%
%
%
%
% Olavo Badaro Marques, 19/Apr/2018.


%%

rho0 = 1025;


%%

KEphasor = u .* conj(v);


%%

%
ctefactor = 0.5 * rho0/2;

%
KE = ctefactor * (UmAmp.^2 + VmAmp.^2);

