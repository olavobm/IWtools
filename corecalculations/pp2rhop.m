function [rhop, zoutrhop] = pp2rhop(zpp, pp, z4rhop)
% [rhop, zoutrhop] = PP2RHOP(zpp, pp, z4rhop)
%
%   inputs
%       - zpp: vector of depths (positive) where pp is given.
%       - pp: hydrostatic pressure perbutation.
%       - z4rhop (optional):
%
%   outputs
%       - rhop: density perturbation.
%
%
%
% Usually, zpp is in meters, pp is in Pascals and rhop in kg m^{-3}.
%
% TO DO:
%   - Include different ways to compute the derivative.
%
% See also: rhop2eta.m, pp2eta.m, eta2pp.m.
%
% Olavo Badaro Marques, 11/Apr/2018.


%%

g = 9.8;


%%

if iscolumn(zpp)
	% that's right -- do nothing.
elseif isrow(zpp)
    zpp = zpp(:);
else
    error('zpp must be a vector of numbers.')
end


%%

if length(zpp)~=size(pp, 1)
    error('Length of zpp is not the same as the number of rows of pp.')
end

%
zpp = repmat(zpp, 1, size(pp, 2));


%%
if exist('z4rhop', 'var')
	linterp = true;
else
    linterp = false;
end


%%
% -------------------------------------------------------------------------
% -------------------------- THE CALCULATIONS -----------------------------
% -------------------------------------------------------------------------


%%

dpdz = (pp(1:end-1, :) - pp(2:end, :)) ./ ...
       (zpp(1:end-1, :) - zpp(2:end, :));

zoutrhop = (zpp(1:end-1, 1) + zpp(2:end, 1)) ./ 2;


%%

rhop = dpdz / g;    % there is no minus sign because
                    % depth increases downward


%%

if linterp
	rhop = interp1overnans(zoutrhop, rhop, z4rhop);
end



