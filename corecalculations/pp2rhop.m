function rhop = pp2rhop(pp, zRhoN)
%
%
%   inputs
%       -
%       -
%
%   outputs
%       - rhop: density perturbation.
%
%
%
% Olavo Badaro Marques, 11/Apr/2018.


%%

g = 9.8;


%%

dpdz = pp;


%%

rhop = - dpdz / g;



