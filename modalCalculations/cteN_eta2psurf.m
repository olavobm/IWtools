function psurf = cteN_eta2psurf(N2, D, nmd, etaAmp)
% psurf = CTEN_ETA2PSURF(N2, D, nmd, etaAmp)
%
%   inputs
%       - N2: buoyancy frequency squared (radians per second squared).
%       - D: water depth (in meters).
%       - nmd: mode number.
%       - etaAmp: isopycnal displacement (in meters).
%
%   outputs
%       - psurf: pressure at the surface.
%
% Compute pressure at the surface correspondent to an isopycnal
% displacement etaAmp for normal modes in constant stratification
% (where the displacement modes are sines -- etaAmp is the amplitude
% of the sines). Positive isopycnal displacements are associated
% with negative surface displacements.
%
% Divide by 10000 (1e4) to (approximately) convert psurf to meters.
%
% The equation is obtained by integrating the hydrostatic equation
% and requiring the depth-integral of pressure to be zero.
%
% Olavo Badaro Marques, 02/Nov/2017.

%
rho0 = 1025;

%
psurf = - rho0 .* N2 .* D .* etaAmp ./ (nmd*pi);