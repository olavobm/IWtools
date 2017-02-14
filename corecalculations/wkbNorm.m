function [zwkb, xwkb, zStretch] = wkbNorm(z, N, zx, x, zgrid, N0)
% [zn, xn] = WKBNORM(z, N, x)
%
%   inputs:
%       - z:
%       - N: buoyancy frequency, in radians per second.
%       - zx: depth where x is specified.
%       - x: vector or matrix. The normalization is done for
%            each column separately.
%       - N0 (optional): reference buyancy frequency
%                        (default is the average).
%
%   outputs:
%       - zn: stretched vertical coordinate.
%       - xn: WKB-normalized x.
%
% WKBNORM does the normalization of ....
%
% Velocity only???
% How small can N really be??
%
% Olavo Badaro Marques, 10/Feb/2017.

%%

if isrow(z)
	z = z';
end


%% Make sure there is no N <= 0:

lbadN = (N(:) <= 0);

if any(lbadN)
    
    N(lbadN) = min(N(~lbadN));
    warning('there are N samller or equal to 0.')
%     error('Given bouyancy frequency has values less or smaller than 0.')
end


%% Make sure the first value is given at the surface

% what if is is a matrix??
if z(1)~=0
    error('First element MUST be 0')
end


%% If reference value is not given, set it as the mean of N:

if ~exist('N0', 'var')
    N0 = mean(N(:));
end


%% Compute the stretched vertical coordinate:

zStretch = wkbStretch(z, N, N0);


%% Do the WKB normalization:

% Interpolate N where values of x are given:
Natx = interp1overnans(z, N, zx);

% Normalization itself:
xNorm = wkbScale(x, Natx, N0);

%
zwkb = interp1overnans(z, zStretch, zx);


%% If zgrid is given, interpolate on grid of stretched depths:

if ~exist('zgrid', 'var')
    
    xwkb = xNorm;
   
else
    
    xwkb = NaN(length(zgrid), size(x, 2));
    
    for i = 1:size(x, 2)
    
        xwkb(:, i) = interp1(zwkb, xNorm(:, i), zgrid);
    
    end
    
end



