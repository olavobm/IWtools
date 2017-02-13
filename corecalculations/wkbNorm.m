function [zn, xn] = wkbNorm(z, N, zx, x, N0)
% [zn, xn] = WKBNORM(z, N, x)
%
%   inputs:
%       - z:
%       - N: buyancy frequency, in radians per second.
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
% What do I really need???
% Velocity only???
% How small can N really be??
% What if I integrate upwards, rather than downwards??
%
% Olavo Badaro Marques, 10/Feb/2017.


%% If x or z are a rows vector, transpose them:

if isrow(z)
    z = z';
end

if isrow(x)
    x = x';
end

if isrow(N)
    N = N';
end


%% Make sure there is no N <= 0:

lbadN = (N(:) <= 0);

if any(lbadN)
    
    N(lbadN) = min(N(~lbadN));
    warning('there are N samller or equal to 0.')
%     error('Given bouyancy frequency has values less or smaller than 0.')
end


%% Make sure the first value is given at the surface

if z(1)~=0
    error('First element MUST be 0')
end


%% If reference value is not given, set it as the mean of N:

if ~exist('N0', 'var')
    N0 = mean(N(:));
end


%% Compute the stretched vertical coordinate:

zn = cumtrapz(z, N/N0);


%% Do the WKB normalization:

% Interpolate N where values of x are given:
Natx = interp1overnans(z, N, zx);

% Normalization itself:
xnorm = x ./ repmat(sqrt((Natx/N0)), 1, size(x, 2));

%
znorm = interp1overnans(z, zn, zx);


%% Now interpolate the normalized variable
% onto the stretched coordinate:

xn = NaN(size(zn, 1), size(x, 2));

for i = 1:size(x, 2)
    
    xn(:, i) = interp1(znorm, xnorm(:, i), zn);
    
end


