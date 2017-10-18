function [xyRay] = raytraceplane(x, y, dydx, xy0, ang0, dxStep, nsteps, xyrto)
% [xyRay] = RAYTRACEPLANE(x, y, dydx, xy0, ang0, dxStep, xyrto)
%
%   inputs
%       - x:
%       - y:
%       - dydx:
%       - xy0:
%       - ang0:
%       - dxStep:
%       - nsteps:
%       - xyrto (optional):
%
%   outputs
%       - xzRay: Nx2 with N coordinates of the ray. The first row is xy0.
%
%
%
%
% Olavo Badaro Marques, 18/Oct/2017.


%%

if not((xy0(1)>=x(1) && xy0(1)<=x(end)) && (xy0(2)>=y(1) && xy0(2)<=y(end)))
    error('Initial point is not inside the domain defined by x and y.')
end


%%

if ~exist('xyrto', 'var')
    xyrto = 1;
end


%%

% % [xg, yg] = meshgrid(x, y);


%%

rayQuad = [cos(ang0)/abs(cos(ang0)), sin(ang0)/abs(sin(ang0))];


%% Pre-allocate space for output variable

%
if ~exist('nsteps', 'var') || isempty(nsteps)
    nsteps = round(sqrt(length(x)^2 + length(y)^2));
end

%
xyRay = NaN(nsteps+1, 2);


%%

xyNow = xy0;
indrow = 1;

xyRay(indrow, :) = xyNow;


%%


for i = 1:nsteps
    
    %% Interpolate dydx to the current point
    
    dydxNow = interp2(x, y, dydx, xyNow(1), xyNow(2));
    
    
    %%

    % dydxNow is positive and the output of atan function (always
    % in the first quadrant). THEREFORE, we must change rayAng to
    % be in the right quadrant.
    rayAng = dydxNow;
    
    % Now change the ray angle (reflect it across cartesian
    % axes) depending on the rayQuad, which determines the
    % quadrant the ray is:
    if rayQuad(1) < 0
        rayVec = reflect2Dacrossline([0; 0], [0; 1], ...
                                     [cos(rayAng); sin(rayAng)]);
                                 
        rayAng = atan2(rayVec(2), rayVec(1));
    end
	
    if rayQuad(2) < 0
        rayVec = reflect2Dacrossline([0; 0], [1; 0], ...
                                     [cos(rayAng); sin(rayAng)]);
                                 
        rayAng = atan2(rayVec(2), rayVec(1));
    end
    
    
    %% Define trace step by scaling the horizontal distance step
    
    traceStep = dxStep;
    
    
    %% Trace next point on the ray:
    xyTrc(1) = xyNow(1) + (traceStep .* cos(rayAng));
    xyTrc(2) = xyNow(2) + (traceStep .* sin(rayAng));
       
    xyNow = xyTrc;
    
    
    %%
    
    xyRay(i+1, :) = xyNow;
    
    
    %% If current point is outside the domain, then break the loop
    
    if not((xyNow(1)>=x(1) && xyNow(1)<=x(end)) && ...
           (xyNow(2)>=y(1) && xyNow(2)<=y(end)))
       
        break
        
    end
    
    
end


%% Remove NaN points added when pre-allocating

lraypts = ~isnan(xyRay(:, 1));

xyRay = xyRay(lraypts, :);


       
    
    



