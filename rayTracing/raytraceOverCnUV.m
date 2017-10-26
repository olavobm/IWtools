function [xyRay, cnRay] = raytraceOverCnUV(lon, lat, cn, U, V, xya0, dtN)
% [xyRay, cnRay] = RAYTRACEOVERCNUV(lon, lat, cn, U, V, xya0, dtN)
%
%   inputs
%       - lon: longitude vector of the domain.
%       - lat: latitude    "    "   "     "
%       - cn: eigenspeed field (for every lon/lat coordinate).
%       - U:
%       - V:
%       - xya0: 1x3 array with initial x/y positions and direction.
%       - dtN: 1x2 array with time resolution and
%              total number of time steps.
%
%   outputs
%       - xzRay: Nx2 with N coordinates of the ray. The first row is xy0.
%       - cnRay: Nx1 array with eigenspeeds along the ray.
%
%
% Olavo Badaro Marques, 25/Oct/2017.


%%

wvfreq = 2*pi / (12.42*3600);

%
dt = dtN(1);
nsteps = dtN(2);

%
xyRay = NaN(nsteps+1, 2);
cnRay = NaN(nsteps+1, 1);


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


%% Compute the spatial derivatives of the eigenspeed field

%
cn_x = NaN(Nlat, Nlon);
cn_y = NaN(Nlat, Nlon);

%
dlon_1 = distance(latg(:, 1:end-2), long(:, 1:end-2), ...
                  latg(:, 3:end), long(:, 3:end));
dlon_2 = distance(latg(:, 1), long(:, 1), ...
                  latg(:, 2), long(:, 2));
dlon_3 = distance(latg(:, end-1), long(:, end-1), ...
                  latg(:, end), long(:, end));
              
dlon_1 = 1000 * deg2km(dlon_1);
dlon_2 = 1000 * deg2km(dlon_2);
dlon_3 = 1000 * deg2km(dlon_3);

cn_x(:, 2:end-1) = (cn(:, 3:end) - cn(:, 1:end-2)) ./ dlon_1;
cn_x(:, 1)   = (cn(:, 2) - cn(:, 1)) ./ dlon_2;
cn_x(:, end) = (cn(:, end) - cn(:, end-1)) ./ dlon_3;

%
dlat_1 = distance(latg(1:end-2, :), long(1:end-2, :), latg(3:end, :), long(3:end, :));
dlat_2 = distance(latg(1, :), long(1, :), latg(2, :), long(2, :));
dlat_3 = distance(latg(end-1, :), long(end-1, :), latg(end, :), long(end, :));
              
dlat_1 = 1000 * deg2km(dlat_1);
dlat_2 = 1000 * deg2km(dlat_2);
dlat_3 = 1000 * deg2km(dlat_3);

cn_y(2:end-1, :) = (cn(3:end, :) - cn(1:end-2, :)) ./ dlat_1;
cn_y(1, :)   = (cn(2, :) - cn(1, :)) ./ dlat_2;
cn_y(end, :) = (cn(end, :) - cn(end-1, :)) ./ dlat_3;


%% Compute the spatial derivatives of the velocity field

%
U_x = NaN(Nlat, Nlon);
U_y = NaN(Nlat, Nlon);
V_x = NaN(Nlat, Nlon);
V_y = NaN(Nlat, Nlon);

% Spatial derivatives of U
U_x(:, 2:end-1) = (U(:, 3:end) - U(:, 1:end-2)) ./ dlon_1;
U_x(:, 1)   = (U(:, 2) - U(:, 1)) ./ dlon_2;
U_x(:, end) = (U(:, end) - U(:, end-1)) ./ dlon_3;

U_y(2:end-1, :) = (U(3:end, :) - U(1:end-2, :)) ./ dlat_1;
U_y(1, :)   = (U(2, :) - U(1, :)) ./ dlat_2;
U_y(end, :) = (U(end, :) - U(end-1, :)) ./ dlat_3;

% Spatial derivatives of V
V_x(:, 2:end-1) = (V(:, 3:end) - V(:, 1:end-2)) ./ dlon_1;
V_x(:, 1)   = (V(:, 2) - V(:, 1)) ./ dlon_2;
V_x(:, end) = (V(:, end) - V(:, end-1)) ./ dlon_3;

V_y(2:end-1, :) = (V(3:end, :) - V(1:end-2, :)) ./ dlat_1;
V_y(1, :)   = (V(2, :) - V(1, :)) ./ dlat_2;
V_y(end, :) = (V(end, :) - V(end-1, :)) ./ dlat_3;


%%

% ------------------------------------------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------


%%

%
xyNow = [xya0(1), xya0(2)];
rayAng = xya0(3);

%
cnpt = interp2(lon, lat, cn, xyNow(1), xyNow(2));
[cppt, cgpt] = cn2cpcg(cnpt, wvfreq * 24*3600/(2*pi), xyNow(2));

%
pxpyNow = [ cos(rayAng)/cppt, ...
            sin(rayAng)/cppt ];

%
xyRay(1, :) = xyNow;
cnRay(1) = cnpt;


%%

traceStep = (cgpt * dt) / 111000;


%%

for i = 1:nsteps
    
    %% --------------------------------------------------------------------
    % Trace next point on the ray:
    xyTrc(1) = xyNow(1) + (traceStep .* cos(rayAng));
    xyTrc(2) = xyNow(2) + (traceStep .* sin(rayAng));
    
    [xyTrc(2), xyTrc(1)] = reckon(xyNow(2), xyNow(1), traceStep, 90 - (180*rayAng/pi));    
     
    %
    xyTrc(1) = wrapPhase([0, 360], xyTrc(1));
    
    %
    xyNow = xyTrc;
    
    % Assign new coordinates to output variable
    xyRay(i+1, :) = xyNow;
    

    %% --------------------------------------------------------------------

    % If current point is outside the domain, then break the loop
    if not((xyNow(1)>=lon(1) && xyNow(1)<=lon(end)) && ...
           (xyNow(2)>=lat(1) && xyNow(2)<=lat(end)))
       
        warning('Ray left the domain')
        break
    end
    
    
    %%
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %     
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    %% --------------------------------------------------------------------
    
    %
    cnpt = interp2(lon, lat, cn, xyNow(1), xyNow(2));
    [cppt, cgpt] = cn2cpcg(cnpt, wvfreq * 24*3600/(2*pi), xyNow(2));
    
    %
    fpt = interp2(lon, lat, f4ray, xyNow(1), xyNow(2));
    
    %
    bpt = interp2(lon, lat, b4ray, xyNow(1), xyNow(2));
    
    %
    Upt = interp2(lon, lat, U, xyNow(1), xyNow(2));
    Vpt = interp2(lon, lat, V, xyNow(1), xyNow(2));
    
    %
    cnRay(i+1) = cnpt;
        
    
    %%
	traceStep = (cgpt * dt) / 111000;
    
    
    %% --------------------------------------------------------------------
    
    % For angles closer to ZONAL
    if abs(tan(rayAng)) <= 2
        
        %
        dcndy = interp2(lon, lat, cn_y, xyNow(1), xyNow(2));
        dUdy = interp2(lon, lat, U_y, xyNow(1), xyNow(2));
        dVdy = interp2(lon, lat, V_y, xyNow(1), xyNow(2));
        
        % Equation (28) in RP 2006
        dHdy = 2*cnpt*(pxpyNow(1)^2 + pxpyNow(2)^2)*dcndy + ...
                        - 2*(Upt*pxpyNow(1)^2 - pxpyNow(1) + Vpt*pxpyNow(1)*pxpyNow(2))*dUdy + ...
                        - 2*(Vpt*pxpyNow(2)^2 - pxpyNow(2) + Upt*pxpyNow(1)*pxpyNow(2))*dVdy + ...
                        (2*fpt/(wvfreq^2))*bpt;
        
        %
        dldx = 2*(cnpt^2 - Upt^2)*pxpyNow(1) + 2*Upt - 2*Upt*Vpt*pxpyNow(2);
        dldx = 1/dldx;
        
        % Equation (27)
        dpydxNow = - dHdy * dldx;
        
        %
        pxpyNow(2) = pxpyNow(2) + ( (111000) * dpydxNow * (traceStep .* cos(rayAng)) );
        
        %
%         pxpyNow(1) = sign(cos(rayAng)) 
        
        pxpyNow(1) = solve4otherP(fpt, wvfreq, cnpt, Upt, Vpt, pxpyNow(2));

        keyboard
        
% %         if ~isreal(pxpyNow)
% %             keyboard
% %         end
        
	% For angles closer to MERIDIONAL
    else
        
        %
        dcndx = interp2(lon, lat, cn_x, xyNow(1), xyNow(2));
        dUdx = interp2(lon, lat, U_x, xyNow(1), xyNow(2));
        dVdx = interp2(lon, lat, V_x, xyNow(1), xyNow(2));
        
        %
        dHdy = 2*cnpt*(pxpyNow(1)^2 + pxpyNow(2)^2)*dcndx + ...
                        - 2*(Upt*pxpyNow(1)^2 - pxpyNow(1) + Vpt*pxpyNow(1)*pxpyNow(2))*dUdx + ...
                        - 2*(Vpt*pxpyNow(2)^2 - pxpyNow(2) + Upt*pxpyNow(1)*pxpyNow(2))*dVdx;
        
        %
        dldy = 2*(cnpt^2 - Vpt^2)*pxpyNow(2) + 2*Vpt - 2*Upt*Vpt*pxpyNow(1);
        dldy = 1/dldy;
        
        %
        dpxdyNow = - dHdy * dldy;
        
        %
        pxpyNow(1) = pxpyNow(1) + ( 111000 * dpxdyNow * (traceStep .* sin(rayAng)) );

        %        
        pxpyNow(2) = solve4otherP(fpt, wvfreq, cnpt, Vpt, Upt, pxpyNow(1));

        
% %         if ~isreal(pxpyNow)
% %             keyboard
% %         end
        
    end
    

    %% --------------------------------------------------------------------
    
    % Update the ray angle
    rayAng = atan2(pxpyNow(2), pxpyNow(1));

end


end

%%

function p2 = solve4otherP(fpt, wvfreq, cnpt, u1, u2, p1)

    %
    a = cnpt^2 - u1^2;

    b = 2 * u1 * (1 - u2*p1);

    c = (cnpt^2 - u2^2)*p1^2 + 2*u2*p1 - 1 + (fpt/wvfreq)^2;
    
    %
    p2_1 = (- b - sqrt(b^2 - 4*a*c)) / 2*a;
    p2_2 = (- b + sqrt(b^2 - 4*a*c)) / 2*a;
    
    %
    
    p2 = p2_1;
end








