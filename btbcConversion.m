function C = btbcConversion(gradD, ubt, vbt, pp, avgpts)
% C = BTBCCONVERSION(gradD, ubt, vbt, pp, avgpts)
%
%   inputs
%       - gradD: bathymetry gradient, with x/y derivatives
%                as real/imaginary parts.
%       - ubt: zonal barotropic velocity.
%       - vbt: meridional   "       "
%       - pp: pressure perturbation.
%       - avgpts: number of points (in time) to compute wave-average.
%
%   outputs
%       - C: conversion, in W/m2.
%
% TO DO:
%   - modify to take complex number as inputs
%	- allow different types of input size (for timeseries, map, etc)
%
% Olavo Badaro Marques, 01/Nov/2017


%%

if ~isequal(size(ubt), size(vbt)) || ~isequal(size(ubt), size(pp))
    error('Input variables do not have the same size.')
end


%%

up = ubt.*pp;
vp = vbt.*pp;


%%

C = (real(gradD) .* up) + (imag(gradD) .* vp);


%%

if exist('avgpts', 'var')
    N = length(C);
    C = obmBinAvg(1:N, C, avgpts);
end

