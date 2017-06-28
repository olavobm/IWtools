function [cp, cg] = cn2cpcg(ce, lat)
% [cp, cg] = CN2CPCG(ce, lat)
%
%   inputs:
%       - cn:
%       - lat:
%
%   outputs:
%       - cp: horizontal phase speed.
%       - cg: horizontal group velocity magnitude.
%
% The eigenspeed is the geometric mean of phase speed and group velocity
% magnitude. Without rotation, a shallow-water gravity wave has
% horizontal.... phase speed and group velocity magnitude are equal.
% With rotation such that
%
%   cp = 
%
%
%   cg = 
%
% Olavo Badaro Marques, 28/Jun/2017.


cgall(:, i) = IWM2.ce' .* ( (2*pi/(12.42*3600))^2 - ...
                                (sw_f(IWM2.lat))^2 ) /  ...
                                ((2*pi/(12.42*3600))^2);