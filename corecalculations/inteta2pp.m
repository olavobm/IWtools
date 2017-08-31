function [pp, psurf] = inteta2pp(z, eta, N2, zN2, rho0)
% [pp, psurf] = INTETA2PP(z, eta, N2, zN2, rho0)
%
%   inputs
%       - z: vector of the data points depths in meters. Depth is greater
%            than 0 and should be specified in ascending order.
%       - eta: vector or matrix of the isopycnal displacement.
%       - N2: background buoyancy frequency squared. If it is a
%             vector, it is applied everywhen in eta.
%       - zN2 (optional): depth where N2 is specified. If not given, code
%                         assumes N2 is given at the same depth as eta.
%       - rho0 (optional): reference potential density (default is 1025).
%
%   outputs
%       - pp: pressure perturbartion.
%       - psurf: baroclinic surface pressure.
%
%
%
% Olavo Badaro Marques, 31/Aug/2017.