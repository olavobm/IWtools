function [cp, cg] = cn2cpcg(cn, freq, lat)
% [cp, cg] = CN2CPCG(cn, freq, lat)
%
%   inputs:
%       - cn: eigenspeed of the n'th vertical mode.
%       - freq: wave frequency (in radians per second).
%       - lat: latitude (in degrees).
%
%   outputs:
%       - cp: horizontal phase speed (in m/s).
%       - cg: horizontal group velocity magnitude (in m/s).
%
% The eigenspeed is the geometric mean of phase speed (cp) and group
% velocity (cg) magnitude. Without rotation, cp and cg have the same
% magnitude (i.e. due to the shallow-water approximation).
%
% Olavo Badaro Marques, 28/Jun/2017.


%% Calculate Coriolis parameter:

omegaEarth = 7.292115e-5;	% Earth's rotation rate

f0 = 2 * omegaEarth * sin(lat);


%% Compute phase speed and group velocity magnitude:

%
cp = cn .* freq/(sqrt(freq.^2 - f0.^2));

%
cg = cn .* sqrt(freq.^2 - f0.^2)/freq;



