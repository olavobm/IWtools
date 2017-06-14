function gmSpec = plotGM76(f0, N0, fgrid, lnewfig)
% gmSpec = PLOTGM76(f0, N0, fgrid)
%
%   inputs:
%       - f0
%       - N0
%       - fgrid (optional):
%       - lnewfig (optional): default is false.
%
%   outputs:
%       - gmSpec:
%
% TO DO: break this function down into calculating / plotting.
%
% Olavo Badaro Marques, 09/Jun/2017.

%%

if ~exist('lnewfig', 'var')
    lnewfig = false;
end


%%

f0 = abs(f0);
flow = f0;


%%

if ~exist('fgrid', 'var') || fgrid
    
    xlimaux = get(gca, 'XLim');
%     fgrid = linspace(flow, xlimaux(2), 5000);
    fgrid = linspace(flow, 0.001, 5000);
else
    
    fgrid_norm = fgrid / ((24*3600)/(2*pi));
    
end


%%

gmSpec = uGM76(f0, N0, fgrid_norm);


%%

if lnewfig
    figure
end

    %
    loglog(gmSpec.freq, gmSpec.spec)