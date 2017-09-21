function depthintPE = modeAmp2PE(pamp, H, N2, modeProf, nmd)
% depthintPE = modeAmp2PE(pamp, H, N2, ntype)
%
%   inputs
%       - pamp: pressur modal amplitude.
%       - H: water depth (in meters).
%       - N2: buoyancy frequency squared (in radians per s^2).
%       - nmd: mode number.
%       - ntype (optional): default is 0.
%
%   outputs
%       - depthintPE: depth0integrated kinetic energy (in J/m^2).
%
%
% case 0 - arbitrary N2
% case 1 - cte N2 and modal shapes are unscaled sines and cosines.
%
% Olavo Badaro Marques, 21/Sep/2017.


%%

%
if ~isempty(modeProf)
    ntype = 0;
else
    ntype = 1;
end

%
rho0 = 1025;


%%

switch ntype

    case 0
        
        %
        zmod = modeProf{1};
        phimod = modeProf{2};
        
        %
        derivMod = NaN(length(zmod), 1);
        
        derivMod(2:end-1) = (phimod(3:end) - phimod(1:end-2)) ./ ...
                              (zmod(3:end) -   zmod(1:end-2));
        derivMod(1) = (phimod(2) - phimod(1)) /  (zmod(2) - zmod(1));
        derivMod(end) = (phimod(end) - phimod(end-1)) /  (zmod(end) - zmod(end-1));
        
        %
        varaux = (1./N2) .* (derivMod).^2;
        pefactor = trapz(zmod, varaux);
        
     
    case 1
        
        pefactor = ((nmd*pi)^2)/(2*H*N2);
        
end


%% Compute depth-integrated KE

depthintPE = (1/(2*rho0)) * pefactor .* (pamp.^2);
