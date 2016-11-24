function [etafcnz, backgheight] = totalVertDisplacement(time, z, sgth)%, meanwndw, sgthlevels)
%
% [etafcnz, backgheight] = TOTALVERTDISPLACEMENT(time, z, sgth)
% 
% (Potential) density given at fixed depths. This must be a MxN matrix:
% rows representing different depths, columns different times.
%
% There can be no gaps. It is your responsability to interpolate over
% the gaps, subset the data, etc... and provide a density field defined
% everywhere in the domain before calling this function.
%
% Data may not have regular time intervals - if not regular, you must
% be very cautious with interpretations, in particular because the
% observed mean may be very different from the true mean.
%
% Olavo Badaro Marques -- 23/Sep/2016


%% Check there are no NaNs in the input. They are NOT allowed:

if any(any(isnan(sgth)))

    error(['Density input has NaNs. Function ' mfilename ' can ' ...
           'not compute displacement. Remove all NaNs.'])
end


%% Number of time-points:
timepts = size(sgth, 2);


%% If non-existent, create the levels of sgth
% for which we want to compute displacement:

if ~exist('sgthlevels', 'var')
    minsgth = min(min(sgth));
    maxsgth = max(max(sgth));
    sgthlevels = linspace(minsgth, maxsgth, length(z));
end


% I BELIEVE, I HAVE TO IMPROVE THIS PART, SUCH THAT WE WON'T HAVE
% NANS AT THE EDGES WHEN TRANSFORMING BACK.

%% Compute isopycnal height. In other words, invert
% the data: go from density being a function of depth
% (and time) to depth being a function of density:

% Pre-allocate space:
zfcnsgth = NaN(size(sgth));

% Loop through time and compute isopycnal height
% of the density values sgthlevels:
for i = 1:timepts
    zfcnsgth(:, i) = interp1(sgth(:, i), z, sgthlevels);
end


%% Compute the "background" isopycnal height (backgheight), where the
% background time scale is given by meanwndw. By background, I mean that
% the "noise" has been averaged out. Here, noise is referred to
% variability at "small" time scales, which can, in principle, be averaged
% out from a time series of length meanwndw. In fact, this noise is the
% signal this function was written for (i.e. the variability from
% the background):

if exist('meanwndw', 'var')
% if ~isempty(meanwndw)
    keyboard
%     % case 1: using window
%     smooth = nan( size(iso) );
%     for i = 1 : length(time)
%         yday0 = yday(i);
%         igd = find( abs(yday-yday0)<=FP.window/2 );
%         smooth(:, i) = nanmean( iso(:, igd), 2);
%     end
    
else 
    
    % Background is the mean of the entire time series:
    backgheight = nanmean(zfcnsgth, 2);
    
    % Replicate to go from vector to matrix
    backgheight = backgheight * ones(1, length(time));
end


%% The vertical displacement of each isopycnal in sgthlevels
% (etafcnsgth) is the height anomaly from the background. The
% is a minus sign in front because depth increases downward
% but we want positive displacements to be upward:

etafcnsgth = - (zfcnsgth - backgheight);


%% Finally, we invert the results back: go from a
% function of density to a function of depth:

% Pre-allocate space:
etafcnz = NaN(size(sgth));

% Loop through time and compute isopycnal height at the z depth levels:
for i = 1:length(time)

    try
        indok = find(~isnan(zfcnsgth(:, i)));
%         etafcnz(:, i) = interp1(zfcnsgth(indok, i), etafcnsgth(indok, i), z);
        
%         % Why not?????
        etafcnz(:, i) = interp1(backgheight(indok, i), etafcnsgth(indok, i), z);
    catch
        keyboard
    end

end

