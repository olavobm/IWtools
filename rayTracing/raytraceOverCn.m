function [xyRay] = raytraceOverCn(lon, lat, cn, xya0)
% [xyRay] = RAYTRACEOVERCN(lon, lat, cn, xya0)
%
%   inputs
%       - lon: longitude vector of the domain.
%       - lat: latitude    "    "   "     "
%       - cn: eigenspeed.
%       - xya0: 1x3 array with initial x/y positions and direction.
%
%   outputs
%       - xzRay: Nx2 with N coordinates of the ray. The first row is xy0.
%
%
% ------------------------------------------------------------------------
%         BE CAREFUL WITH THE DEFINITION OF THE ANGLE THETA USED!!!!
% ------------------------------------------------------------------------
%
% Olavo Badaro Marques, 18/Oct/2017.


%%

Nlat = length(lat);
Nlon = length(lon);

%
[long, latg] = meshgrid(lon, lat);


%% Compute Coriolis parameter and its derivative (i.e. beta)

%
fvec = gsw_f(lat);

% Beta
omegaEarth = 7.292115e-5;
radiusEarth = 6.4e6;
betavec = 2 * omegaEarth * cosd(lat) ./ radiusEarth;


%%

f4ray = repmat(fvec(:), 1, length(lon));
b4ray = repmat(betavec(:), 1, length(lon));


%% Compute the derivative of the eigenspeed

%
cn_x = NaN(Nlat, Nlon);
cn_y = NaN(Nlat, Nlon);

%
cn_x(:, 2:end-1) = (cn(:, 3:end) - cn(:, 1:end-2)) ./ (long(:, 3:end) - long(:, 1:end-2));

cn_x(:, 1)   = (cn(:, 2) - cn(:, 1)) ./ (long(:, 2) - long(:, 1));
cn_x(:, end) = (cn(:, end) - cn(:, end-1)) ./ (long(:, end) - long(:, end-1));

%
cn_y(2:end-1, :) = (cn(3:end, :) - cn(1:end-2, :)) ./ (latg(3:end, :) - latg(1:end-2, :));

cn_y(1, :)   = (cn(2, :) - cn(1, :)) ./ (latg(2, :) - latg(1, :));
cn_y(end, :) = (cn(end, :) - cn(end-1, :)) ./ (latg(end, :) - latg(end-1, :));

% ------------------------------------------------------------
%       CHECK THE SIGNS OF THE DERIVATIVES!!!!!
% ------------------------------------------------------------


%%





%%

% ------------------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------





