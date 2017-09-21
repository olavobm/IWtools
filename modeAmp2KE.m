function depthintKE = modeAmp2KE(uamp, vamp, H, ntype)
% depthintKE = MODEAMP2KE(uamp, vamp, H, ltype)
%
%   inputs
%       - uamp:
%       - vamp:
%       - H:
%       - ntype (optional): default is 0.
%
%   outputs
%       - depthintKE: depth0integrated kinetic energy (in J/m^2).
%
%
% case 0 - standard normalization (depth integral
%          squared is Kronecker delta)
% case 1 - cte N2 and modal shapes are unscaled sines and cosines.
%
%
% Olavo Badaro Marques, 21/Sep/2017.


%%

if ~exist('ntype', 'var')
	ntype = 0;
end

%
rho0 = 1025;


%%

switch ntype
     
    case 0
        
        kefactor = 1;
        
    case 1
        
        kefactor = H/2;
        
% % %     case 2
% % %         
% % %         kefactor = trapz(zmod, ModShape.^2);
end


%% Compute depth-integrated KE

depthintKE = (rho0/2) * kefactor .* (uamp.^2 + vamp.^2);
