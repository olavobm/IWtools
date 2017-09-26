function depthintKE = modeAmp2KE(uamp, vamp, H, ntype)
% depthintKE = MODEAMP2KE(uamp, vamp, H, ntype)
%
%   inputs
%       - uamp: u modal amplitude amplitude.
%       - vamp: v   "       "         "
%       - H: water depth (in meters).
%       - ntype (optional): normalization type of the modal
%                           shapes (default is 0).
%
%   outputs
%       - depthintKE: time-averaged, depth-integrated
%                     kinetic energy (in J/m^2).
%
% Compute time-averaged, depth-integrated Kinetic Energy from the modal
% amplitude of velocity components. The time-average value is NOT a 
% function of the phase difference between the components.
%
% When taking the depth integral, there is a constant factor that
% multiplies the velocity squared that comes from the integral
% of the modal shapes. The type of normalization is specified by
% the optional input "ntype". The avalable possibilities are:
% 
%       * ntype == 0 - standard normalization (depth integral
%                      of modal shaped squared is one).
%       * ntype == 1 - cte N2 and modal shapes are unscaled cosines.
%
% Olavo Badaro Marques, 21/Sep/2017.


%% Check inputs, set constant

if ~exist('ntype', 'var')
	ntype = 0;
end

%
rho0 = 1025;


%% Define constant factor based on the type of normalization

switch ntype
     
    case 0
        
        kefactor = 1;
        
    case 1
        
        kefactor = H/2;
        
end


%% Compute depth-integrated, time-averaged KE

depthintKE = (1/2) * (rho0/2) * kefactor .* (uamp.^2 + vamp.^2);

