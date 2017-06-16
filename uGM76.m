function specGM = uGM76(f0, N0, fgrid, freqFac)
% specGM = UGM76(f0, N0, fgrid)
%
%   inputs:
%       - f0: Coriolis parameter (radians per second).
%       - N0: bouyancy frequency (radians per second).
%       - fgrid: frequencies (radians per second) where the GM spectrum
%                is computed.
%       - freqFac (optional): ....
%
%   outputs:
%       - specGM: structure with two fields (freq and spec), with the
%                 GM spectrum -- the frequencies and the values of the
%                 spectrum.
%
% Olavo Badaro, 12/Jun/2017.


%%

f0 = abs(f0);


%%

params = Gm76Params;

S = GmOm('Vel', fgrid, f0, N0, params);


%%

if ~exist('freqFac', 'var')
    
    freqFac = 1;
    
else
    
    if ischar(freqFac)
    
        if strcmp(freqFac, 'cpd')
            freqFac = (24*3600)/(2*pi);
        end
    
    else


    end
    
end



%% Convert frequency vector from radians
% per second to cycles per day:

specGM.freq = fgrid * freqFac;
specGM.spec = S / freqFac;