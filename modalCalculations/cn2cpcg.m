function [cp, cg] = cn2cpcg(cn, freq, lat)
% [cp, cg] = CN2CPCG(cn, freq, lat)
%
%   inputs:
%       - cn: eigenspeed of the n'th vertical mode.
%       - freq: wave frequency (in cycles per day).
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


%% Calculate Coriolis parameter

omegaEarth = 7.292115e-5;	% Earth's rotation rate

f0 = 2 * omegaEarth * sind(lat);


%% Convert wave frequency from cycles per day to radians per second

freq = freq .* (2*pi/(24*3600));


%% Compute phase speed and group velocity magnitude

times_factor = freq./(sqrt(freq.^2 - f0.^2));

%
cp = cn .* times_factor;

%
cg = cn ./ times_factor;



