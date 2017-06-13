function specGM = uGM76(f0, N0, fgrid)
%
%
%
% Olavo Badaro, 12/Jun/2017.


%%

f0 = abs(f0);


%%

params = Gm76Params;

S = GmOm('Vel', fgrid, f0, N0, params);


%% Convert frequency vector from radians
% per second to cycles per day:

freqFac = (24*3600)/(2*pi);

specGM.fvec = fgrid * freqFac;
specGM.spec = S / freqFac;