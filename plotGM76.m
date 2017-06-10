function gmSpec = plotGM76(f0, N0, fgrid)
% gmSpec = PLOTGM76(f0, N0, fgrid)
%
%   inputs:
%       - f0
%       - N0
%       - fgrid (optional):
%
%   outputs:
%       - gmSpec:
%
% TO DO: break this function down into calculating / plotting.
%
% Olavo Badaro Marques, 09/Jun/2017.


%%

f0 = abs(f0);
flow = f0;


%%

if ~exist('fgrid', 'var')
    
    xlimaux = get(gca, 'XLim');
%     fgrid = linspace(flow, xlimaux(2), 5000);
    fgrid = linspace(flow, 0.001, 5000);
    
end


%%

params = Gm76Params;

S = GmOm('Vel', fgrid, f0, N0, params);


%% Convert frequency vector from radians
% per second to cycles per day:

freqFac = (24*3600)/(2*pi);

fvecGM = fgrid * freqFac;
SpecGM = S / freqFac;


%%

gmSpec.freq = fvecGM;
gmSpec.psd = SpecGM;


%
% figure
    loglog(gmSpec.freq, gmSpec.psd)