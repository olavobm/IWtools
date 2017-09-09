function [xzRay] = raytraceIW(xg, zg, N, f0, wvf, xz0, rayQuad, traceDx, botstruct)
% [xzr] = RAYTRACEIW(xg, zg, N, f0, wvf, xz0, rayQuad, traceDx, botstruct)
%
%   inputs
%       - xg: horizontal grid points.
%       - zg: vertical grid points.
%       - N: matrix of size (length(zg))x(length(xg)) with the buoyancy
%            frequency, in radians per second. NaNs in the matrix are
%            interpreted as solid boundary.
%       - f0: Coriolis parameter in radians per second.
%       - wvf: wave frequency, in radians per second.
%       - xz0: 1x2 vector with initial (x, z) coordinates of the wave group.
%       - rayQuad: 1x2 vector. Its values can either be -1 or +1. The four
%                  combinations indicate which trignometric quadrant.
%       - botstruct: structure with two fields:
%               * x: horizontal location.
%               * z: bottom depth.
%
%   outputs
%       - xzRay: Nx2 with N coordinates of the ray.
%               xzr(1, :) is always equal to xz0.
%
%
% The tangent (dz/dx) of the ray is giving by the dispersion relationship:
% dz/dx = +- sqrt((wvf^2 - f0^2) / (N^2 - wvf^2))
%
% GOT TO FIGURE OUT HOW TO PROPERLY IDENTIFY THE BOTTOM!!!
%
% TO DO:
%   - The grid resolution gives me angle resolution of tracing. Use
%     that as a diagnostic.
%   - Trace "in time"????
%   - Later take a look at asymptotic expansion for turning depths.
%   - It could be nice to deal with wavenumber of fake magnitude (but
%     correct), which could allow me use some linear algebra and write a
%     better code.
%   - MAKE SURE THAT I WON`T HAVE A PROBLEM IF THE RAY GETS TRACED THROUGH
%     THE BOTTOM BUT OUTSIDE OF THE GRID.
%   - I PROBABLY WANT TO INTEGRATE DZ/DX OR DX/DZ DEPENDING
%     ON THE ANGLE OF THE RAY
%   - PROBABLY BE USEFUL TO ALLOW GIVING THE BOTTOM AS INPUT
%
%   - WOULD BE NICE IF THE GRID DID NOT HAVE TO GO ALL THE
%     WAY TO THE SURFACE
%
%   - IF THE POINT IS CLOSE TO THE BOUNDARY, BUT PAST THE LAST USEFUL
%     GRID POINT, I GET AN ERROR BECAUSE N2 IS NAN BUT RAY HASN`T
%     CROSSED THE BOUNDARY. THIS SUGGEST N2 MUST BE SPECIFIED PAST THE
%     BOTTOM AND THEN MY CRITERION FOR BOTTOM ENCOUNTERING SHOULD CHANGE
%
%   - ACTUALLY, MAYBE THIS FUNCTION IS USELESS FOR THIS PURPOSE. BECAUSE IF
%     I CAN DO WKB STRETCHING, THEN I CAN COME UP WITH A MUCH MORE
%     EFFICIENT ALGORITHM. ANYWAY, I COULD USE THIS CODE FOR A GENERAL RAY
%     TRACING ALGORITHM.
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
% SHOULD ADD NANS AS THE FIRST ROW OF N SUCH THAT THE SURFACE CAN BE
% TREATED AS FLAT BOUNDARY??? THE CODE FOR DOING REFLECTION OFF OF THE
% SURFACE AND THE BOTTOM WOULD THEN BE THE SAME


%% Create a Nx2 matrix with all the N grid point coordinates:

[xgmesh, ygmesh] = meshgrid(xg, zg);

gridPoints = [xgmesh(:), ygmesh(:)];


%% Throw error if user did not specify which
% quadrant the ray is propagating in:
% (ADD ERROR ONLY WHEN THE FUNCTION IS DONE)

if ~exist('rayQuad', 'var')
    rayQuad = [1 , 1];
%     error('bla')  % when the function is done, it should give an error.
end


%%


if ~exist('traceDx', 'var')
    
    % traceDx = median(diff(xg));   % use xg, because this is much larger in
                                    % the ocean than zg, and rays propagate
                                    % mostly horizontally
    traceDx = 1000;

end


%% Identify the bottom (still need generalize the code):


if find(isnan(N), 1)  % I'm not sure if the if is necessary.
    domainBottom = defineNaNregions(isnan(N));
    
    %
    for i = 1:length(domainBottom)
        
        bla = false(size(N));
        
        indbot = sub2ind(size(bla), domainBottom{i}(:, 1), domainBottom{i}(:, 2));
        
        bla(indbot) = true;
%         bla(domainBottom{i}(:, 1), domainBottom{i}(:, 2)) = true;
        
        indbdry = findBoundary(bla);
        indbdry2 = sortHoop(indbdry);
    end
    
    [a, b] = meshgrid(xg, zg);
    
    indhoop = sub2ind(size(a), indbdry2(:, 1), indbdry2(:, 2));

    indhoop = {indhoop};
    
end



%%

xzNow = xz0;

indRay = zeros(1, size(gridPoints, 1));
xzRay = NaN(length(indRay), 2);

% Make sure the first point is inside the grid:
[linX, linZ] = checklGrid(xg([1 end]), zg([1 end]), xzNow);

if linX && linZ
    linGrid = true;
else
    linGrid = false;
end


%% Now make sure it is a non-NaN place (IT WOULD ACTUALLY BE GOOD TO ALLOW
% THE FIRST POINT RIGHT ON THE BOUNDARY, A SIMPLE AND FAIRLY GOOD SOLUTION
% IS TO KEEP THE RAY TRACING IF THE POINT JUST ABOVE HAS A VALID N):

% Interpolate a vertical profile
auxNprof = interp1overnans(xg, N', xzNow(1), xg(2)-xg(1));  
auxNprof = auxNprof';

% Interpolate at a depth in the vertical profile:
nowN2 = interp1(zg, auxNprof, xzNow(2));

if isnan(nowN2)
    linGrid = false;
end


%% Loop for the ray tracing:

%
xzTrc = NaN(1, 2);    % this is useles...

%
i = 1;


while linGrid
    %% Check whether current point is on the grid and
    % whether buoyancy frequency is in a non-NaN
    
	[linX, linZ] = checklGrid(xg([1 end]), zg([1 end]), xzNow);
           
    if ~linX || xzNow(2)>zg(end)
        
        % do nothing because it has left the side boundaries
        %
        % Would probably be good to use the interpolation to get
        % the last point on the edge of the grid.
        
        linGrid = false;
%         break
    end
    
    
    %% Interpolate N to the current point:   
    if linX && linZ
        
        
        % Interpolate a vertical profile
        auxNprof = interp1overnans(xg, N', xzNow(1), xg(2)-xg(1));  
        auxNprof = auxNprof';
        
        % Interpolate at a depth in the vertical profile:
        nowN2 = interp1(zg, auxNprof, xzNow(2));
        
    end
    
    
    %% Reflection off the bottom or the surface:
    % (I PROBABLY WANT TO SPLIT THIS, DEALING WITH POINTS OUTSIDE
    % OF THE GRID FIRST, AND BOTTOM REFLECTION LATER)
    if (isnan(nowN2) && linX && linZ) || (xzNow(2)<zg(1) && linX)
        
        if isnan(nowN2)
            
            % (xzNow(1) - xzRay(i-1, 1))
            
            % This finds the bottom:
            llookforbot = true;
            indlookbot = 1;
            while llookforbot
                
                % there is a problem if the previous point is on (very
                % close to the boundary). For now, I'll use a different
                % point on the ray (i-2) rather than (i-1)).
                %
                % STEPPY BOUNDARY MAY GIVE MULTIPLE INTERSECTIONS --
                % THIS FORCES ME TO USE A SMOOTHED LOW-RES BOUNDARY
                try
                [xcross, ycross, indAtHoop] = intersections(a(indhoop{indlookbot}), ...
                                                            b(indhoop{indlookbot}), ...
                                                            [xzNow(1), xzRay(i-2, 1)], ...
                                                            [xzNow(2), xzRay(i-2, 2)]);
                catch
                    keyboard
                end
                % --------------------------------------------------
                % REMOVE THIS WHEN I ADDRESS NEAR-BOTTOM
                % BUOYANCY FREQUENCY!!!!!
                if isempty(xcross)
                    dxaux = xzNow(1) - xzRay(i-2, 1);
                    dzaux = xzNow(2) - xzRay(i-2, 2);
                    
                    xzNow(1) = xzNow(1) + dxaux;
                    xzNow(2) = xzNow(2) + dzaux;
                    
                    [xcross, ycross, indAtHoop] = intersections(a(indhoop{indlookbot}), ...
                                                            b(indhoop{indlookbot}), ...
                                                            [xzNow(1), xzRay(i-2, 1)], ...
                                                            [xzNow(2), xzRay(i-2, 2)]);
                    
                end
                % --------------------------------------------------
                
                if isempty(xcross)
                    
                    indlookbot = indlookbot + 1;
                    
                else
                    
                    xcross = xcross(1);
                    ycross = ycross(1);  % remove this in the future
                    
                    botstruct.x = a(indhoop{indlookbot});
                    botstruct.z = b(indhoop{indlookbot});
                    
                    llookforbot = false;
                end 
            end
            
            xaux = xcross;
            yaux = ycross;
            
            % If bottom is given I may or may not want to smooth the
            % bottom. If not given, I definitely want to smooth it.
            % THIS CAN GIVE HUGE PROBLEMS AT THE EDGES.... I PROBABLY JUST
            % I WANT TO SPECIFY THE BOTTOM AS INPUT...
            lsmoothbot = true;
            
            if lsmoothbot
                
                smoothscale = 5000;   % 5 km
                
                upperind = ceil(indAtHoop);
                lowerind = floor(indAtHoop);
                
                lwidenpts = true;
                while lwidenpts
                    
                    distSmooth = sqrt((botstruct.x(upperind) - botstruct.x(lowerind)).^2 + ...
                                      (botstruct.z(upperind) - botstruct.z(lowerind)).^2);
                                  
                    if distSmooth < smoothscale
                        lowerind = lowerind - 1;
                        upperind = upperind + 1;
                    else
                        lwidenpts = false;
                    end           
                end
                
                xslope = botstruct.x([lowerind upperind]);
                zslope = botstruct.z([lowerind upperind]);
                
                if xslope(1)>xslope(2)
                    inddslope = [1, 2];
                elseif xslope(1)>xslope(2)
                    inddslope = [2, 1];
                else
                    inddslope = [1, 2];
                end
                
                dxslope = xslope(inddslope(1)) - xslope(inddslope(2));
                dzslope = zslope(inddslope(1)) - zslope(inddslope(2));
                
                dslope = dzslope / dxslope;
                % this should be outside the if block, because I
                % also want to do if user specifies bathymetry
                
                % this is only half of the trigonometric circle.
                % SHOULD I MAKE IT ALL AROUND??? PROBABLY NOT NECESSARY
                
            end
            
            
            % if ray is steeper than slope -- reflect vertical
            % direction
            %
            % if slope is steeper than ray -- reflect horizontal
            % direction
            if abs(tan(rayAng)) > abs(dslope)
                rayQuad(2) = - rayQuad(2);
            else
                rayQuad(1) = - rayQuad(1);
            end
       
%             % Distance between previous and current points:
%             rayPtsDist = sqrt((xzNow(1) - xzRay(i-1, 1))^2 + ...
%                               (xzNow(2) - xzRay(i-1, 2))^2);
%                           
%             % Look in the neighbouring region (i.e. within rayPtsDist)
%             % for the boundary .... MAYBE I COULD IMPLEMENT THIS LATER
             

        else
            % Reflection off the surface:
            
            [xaux, yaux] = intersections([xzRay(i-1, 1), xzNow(1)], ...
                                         [xzRay(i-1, 2), xzNow(2)], ...
                                         [xg(1), xg(end)], [zg(1), zg(1)]);
            
            % Direct the ray towards the ocean (increasing depth):
            rayQuad(2) = 1;
                                     
        end
        
%         if isempty(xaux)
%             keyboard
%         end
        
        xzNow = [xaux, yaux];
        
        % Interpolate N at the current location:
        auxNprof = interp1overnans(xg, N', xzNow(1), xg(2)-xg(1));  
        auxNprof = auxNprof';
        
        % Interpolate at a depth in the vertical profile:
%         nowN2 = interp1(zg, auxNprof, xzNow(2));
        if isnan(nowN2)
            auxbla = auxNprof(~isnan(auxNprof));
            nowN2 = auxbla(end);  % approximate with the
                                  % closest to the bottom
        else
            nowN2 = auxNprof(1);   % surface value
        end
        
        if isnan(nowN2)
            error('nowN2 is NaN. this should not happen here!')
        end
        
    end
    
    
    %% Assign current point to output:
    
	xzRay(i, :) = xzNow;

    i = i + 1;
    
    
    %% Now use N at xzNow to trace to the next point:
    
    % Compute first-quadrant tangent of the internal wave characteristic:
    rayTan = abs( sqrt((wvf^2 - f0^2)/(nowN2^2 - wvf^2)) );
    
    rayAng_1stQuad = atan(rayTan);  % the above with always gives an
                                    % angle in the first quadrant
    rayAng = rayAng_1stQuad;
    
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
       
    
    %% Trace next point on the ray:
    xzTrc(1) = xzNow(1) + traceDx .* cos(rayAng);
    xzTrc(2) = xzNow(2) + traceDx .* sin(rayAng);
       
    xzNow = xzTrc;
          	
    
%     if linX && linZ
%        
% %         indRay(i) = indclosest;   % probably useless now
%         
%         xzRay(i, :) = xzTrc;
% 
%         xzNow = xzTrc;
%         i = i + 1;
%         
%         
%     else
%         
%         % For now, break tracing when ray leaves the
%         % domain THROUGH ANY SIDES OF DOMAIN
%         
%         break
%         
%         
%     end
           
    
end


%%

xzRay = xzRay(~isnan(xzRay(:, 1)), :);



end


%% -----------------------------------------------------------------
% ------------------------------------------------------------------
% ------------------------------------------------------------------

function [linX, linZ] = checklGrid(xglims, zglims, xzpt)
    %%
    
    linX = (xzpt(1) >= xglims(1)) && (xzpt(1) <= xglims(2));
    
    linZ = (xzpt(2) >= zglims(1)) && (xzpt(2) <= zglims(2));
    


end


% function interpNatxz()
%     %% TOTALLY MAKES MY LIFE EASIER TO WRITE THIS!
% 
%     
%     
% end

