function etabt = bt2eta(U, s, D, z, wvfreq, ti)
% etabt = BT2ETA(U, s, D, z, wvfreq, ti)
%
%   inputs
%       - U: upslope barotropic velocity (parallel
%            to the bathymetry gradient).
%       - s: (tangent) of the bottom slope.
%       - D: water depth.
%       - z: depth points where the output is computed.
%       - wvfreq: wave frequency (in radians per second).
%       - ti: time correspondent to U.
%
%   outputs
%       - etabt: isopycnal displacement due to barotropic
%                flow over a 1D linear slope.
%
% BT2ETA computes the isopycnal displacement caused by a tidal
% barotropic flow across linearly sloping topography.
%
%   TO DO:
%       - Check calculcation with complex numbers.
%
% Olavo Badaro Marques, 06/Nov/2017.


%%

if exist('ti', 'var')

    coefaux = - (diff(U) ./ diff(ti)) ./ (wvfreq^2);
    
else
   
    coefaux = U/wvfreq * exp(-1i.*pi/2);
    
end


% Compute the isopycnal displacement
etabt =  s .* (z(:)./D) * coefaux;
