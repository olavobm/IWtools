function [varout] = modeflux2var(F, D, N2, lat0, freq, nmd, listvar)
% [varout] = MODEFLUX2VAR(F, D, N2, lat0, freq, nmd, listvar)
%
%   inputs
%       - F: depth-integrated energy flux (in W/m).
%       - D: water depth.
%       - N2: buoyancy frequency squared (radians per second squared).
%       - lat0: latitude.
%       - freq: wave frequency (cycles per day).
%       - nmd (optional): mode number (default is 1).
%       - listvar (optional): cell array with string of
%                             variables to look at (default is u).
%
%   outputs
%       - varout:
%
%
% Result for constant buoyancy frequency.
%
% ????? AM I SURE OF THESE EQUATIONS ?????
%
%   TO DO:
%       - clean the code.
%       - add comments.
%
% Olavo Badaro Marques, 20/Sep/2017.


%%

if ~exist('nmd', 'var') || isempty(nmd)
	nmd = 1;
end

if ~exist('listvar', 'var')
	listvar = {'u'};
else
    if ~iscell(listvar)
        listvar = {listvar};
    end
end

%
rho0 = 1025;


%%

f0 = gsw_f(lat0);


%%

freq = ((2*pi)/(24*3600)) * freq;


%% 

for i1 = 1:length(nmd)
    
    %
    kmag = (nmd*(pi/D)) * sqrt( (freq^2 - f0^2)/(N2 - freq^2) );
    
    %
    for i2 = 1:length(listvar)
        
        switch listvar{i2}
            
            case 'u'
                
                varout = F * 4 * rho0 * (freq^2 - f0^2) / (freq * kmag * D);
                
                varout = sqrt(varout);
                
                varout = varout * ((freq*kmag) / (rho0 * (freq^2 - f0^2)));
                
            case 'p'
                
                varout = F * 4 * rho0 * (freq^2 - f0^2) / (freq * kmag * D);
                
                varout = sqrt(varout);
                
        end
        
    end
end

