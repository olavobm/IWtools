function [xr, zr] = raytraceIW(xg, zg, N, f0, wvf, xz0, rayQuad)
% [xr, zr] = RAYTRACEIW(xg, zg, N2, f0, wvf, xz0)
%
%   inputs:
%       - xg: horizontal grid points.
%       - zg: vertical grid points.
%       - N: matrix of size (length(zg))x(length(xg)) with the buoyancy
%            frequency, in radians per second. NaNs in the matrix are
%            interpreted as solid boundary.
%       - f0: Coriolis parameter in radians per second.
%       - wvf: wave frequency, in radians per second.
%       - xz0: 1x2 vector with initial (x, z) coordinates of the wave group.
%       - pointDir:
%
%   outputs:
%       - xr: horizontal points along the ray.
%       - zr: vertical points along the ray.
%
%
%
% The tangent (dz/dx) of the ray is giving by the dispersion relationship:
% dz/dx = +- sqrt((wvf^2 - f0^2) / (N^2 - wvf^2))
%
% TO DO:
%   - The grid resolution gives me angle resolution of tracing. Use
%     that as a diagnostic.
%   - Trace "in time"????
%   - Later take a look at asymptotic expansion for turning depths.
%   - It could be nice to deal with wavenumber of fake magnitude (but
%     correct), which could allow me use some linear algebra and write a
%     better code.
%   - DISTANCE, X AND Z WITH DIFFERENT SCALES (!!!!)
%
% The while loop may make the code more synthetic if I first store the ray
% location and then trace it, leaving the next point for subsequent
% iteration.
%
% Olavo Badaro Marques, 14/Feb/2017.

% There is some arbitrariness on the ray location with respect to grid
% points. Some possibilities are:
%   - Trace in time and interpolate N at every time step. That may be good
%  	choosing a small time step, but a long one will give wrong results.
%   - Approximate ray locations by the grid points (what I'm currently
%   doing, but then I have to case about how far the rays are from grid
%   points).
%
% TRACING "IN TIME" IS THE RIGHT APPROACH!!!, BUT THINK OF IT AS DISTANCE!


%% Create a Nx2 matrix with all the N grid point coordinates:

[xgmesh, ygmesh] = meshgrid(xg, yg);

gridPoints = [xgmesh(:), ygmesh(:)];


%%

traceDx = max([median(diff(xg)), median(diff(zg))]);

% allow this to be given as inpute


%%

xzNow = xz0;

linGrid = (xzNow(1)>=xg(1) && xzNow(1)>=xg(end) && ...
           xzNow(2)>=zg(1) && xzNow(2)>=zg(end));

indRay = zeros(1, size(gridPoints, 1));

indclosest = dsearchn(gridPoints, xzNow(:));
    
indRay(1) = indclosest;
    
xzNow = gridPoints(indclosest, :);

xzTrc = NaN(1, 2);
%
i = 2;

while linGrid
    
    rayTan = abs( sqrt((wvf^2 - f0^2)/(N(indclosest)^2 - wvf^2)) );
    
    rayAng = atan(rayTan);
    
%     if rayQuad(1)<0
%         
%     end
    xzTrc(1) = xzNow(1) + traceDx .* cos(rayAng);
    xzTrc(2) = xzNow(2) + traceDx .* sin(rayAng);
    
    indclosest = dsearchn(gridPoints, xzNow(:));
    
    linGrid = (xzNow(1)>=xg(1) && xzNow(1)>=xg(end) && ...
               xzNow(2)>=zg(1) && xzNow(2)>=zg(end));
    
    
    indRay(i) = indclosest;
    
    
    
    i = i + 1;
end










