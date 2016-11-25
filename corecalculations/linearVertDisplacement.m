function [eta, backgrnd] = linearVertDisplacement(z, x, cutoff, wnd, t)
% [eta, backgrnd] = LINEARVERTDISPLACEMENT(z, x, meanwndw, t)
%
%   inputs:
%       - z: vector (would be nice to implement matrix)
%       - x: REGULARLY SPACED IN TIME (would be nice to relax that, which
%            is basically relaxing the running mean function).
%       - cutoff: returns NaN/0 if dTdz is too small.
%       - wnd (optional): depending on the version of obmRunMean, wnd must
%                         be odd.
%       - t (optional): (after I implement improvement in obmRunMean.
%
%   outputs:
%       - eta:
%       - backgrnd:
%
% Function LINEARVERTDISPLACEMENT computes the vertical (across rows)
% displacement under a linear approximation. This displacement is defined
% by how much it is required to vertically advect the background field,
% with locally constant vertical gradient, in order to explain the
% anomalies of x relative to the background.
%
% What do I do for data of station-like format???
%
% Olavo Badaro Marques, 23/Nov/2016.


%% Check inputs:

% Make sure x is a matrix:
if isvector(x)    
    error('Does not make sense!')
end


%% Apply running mean to compute background x field (transpose
% such that filtering is applied along the rows of x):

backgrnd = obmRunMean(x', wnd);
backgrnd = backgrnd';


%% Compute vertical gradient of x (looping
% through columns to deal with NaNs):

dxdz = NaN(size(x));

for i = 1:size(x, 2)
    
    indok = find(~isnan(x(:, i)));
    dxdz(indok, i) = centeredDeriv(z(indok), x(indok, i));
end


%% Apply running mean to dxdz:

dxdz = obmRunMean(dxdz', wnd);
dxdz = dxdz';


%% Compute vertical displacement:

eta = (x - backgrnd) ./ dxdz;

