function C = btbcConvert(gradD, ubt, vbt, pp)
% C = BTBCCONVERT(gradD, ubt, vbt, pp)
%
%   inputs
%       - gradD:
%       - ubt:
%       - vbt:
%       - pp:
%
%   outputs
%       - C: conversion
%
% TO DO:
%   - modify to take complex number as inputs
%
% Olavo Badaro Marques, 01/Nov/2017


%%

if ~isequal(size(ubt), size(vbt)) || ~isequal(size(ubt), size(pp))
    error('Input variables do not have the same size.')
end


%%

upavg = ubt.*pp;
vpavg = vbt.*pp;


%%

C = (real(gradD) .* upavg) + (imag(gradD) .* vpavg);

