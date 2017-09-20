function [varout] = modeflux2var(F, D, N2, freq, lat0, mdn, listvar)
% [varout] = MODEFLUX2VAR(F, D, N2, freq, lat0, nmd, listvar)
%
%   inputs
%       - F: depth-integrated energy flux (in W/m).
%       - D: water depth.
%       - N2: buoyancy frequency squared (radians per second squared).
%       - freq: wave frequency (cycles per day).
%       - lat0: latitude.
%       - mdn (optional): mode number (default is 1).
%       - listvar (optional): cell array with string of
%                             variables to look at (default is u).
%
%   outputs
%       - varout:
%
%
% Result for constant buoyancy frequency.
%
% Olavo Badaro Marques, 20/Sep/2017.


%%

if ~exist('mdn', 'var') || isempty(mdn)
	mdn = 1;
end

if ~exist('listvar', 'var')
	listvar = {'u'};
end


%% 

for i1 = 1:length(mdn)
    
    for i2 = 1:length(listvar)
        
        switch listvar{i2}
            
            case 'u'
                
                
                
            case 'p'
                
                
        end
        
    end
end

