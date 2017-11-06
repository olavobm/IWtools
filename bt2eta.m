function etabt = bt2eta(U, s, D, z, tfreq, ti)
% etabt = BT2ETA(U, s, D, z)
%
%   inputs
%       - U: barotropic velocity perpendicular to the bathymetry.
%       - s: (tangent) of the bottom slope.
%       - D: water depth.
%       - z: depth points to compute the calculation
%       - tfreq:
%       - ti:
%
%   outputs
%       - etabt: isopycnal displacement due to barotropic
%                flow over a 1D linear slope.
%
%
% Olavo Badaro Marques, 06/Nov/2017.


%%

if exist('ti', 'var')

    coefaux = - (diff(U) ./ diff(ti)) ./ (tfreq^2);
    
else
   
    coefaux = U/tfreq;
    
end


%
etabt =  s .* (z(:)./D) * coefaux;
