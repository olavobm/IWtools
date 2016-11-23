function [davg, xavg] = waveAvgLSqrs(x, d, prd, wnd, wndptslide)
% [davg, xavg] = WAVEAVGLSQRS(x, d, prd, wnd, wndptslide)
% 
%   inputs:
%       - x: vector or matrix.
%       - d: variable to be wave-averaged (each of its rows).
%       - prd: period of the sinusoidal, in the same units as x.
%       - wnd (optional): window to, in the same units as x.
%       - wndptslide (optional): this can only be input if wnd is.
%
%   outputs:
%       - davg: wave-averaged variable.
%       - xavg:
%
% Use least squares to compute the wave-average (davg) of d. The
% coefficients of a sinusoidal of period prd plus a constant are
% estimated from least squares and the constant is the wave-averaged
% output is the constant.
% The fit is computed by the function myleastsqrs.m.
%
% If no window is specified, the output is just a number or vector with
% the wave-averaged quantities.
%
% Olavo Badaro Marques, 17/08/2016.


%% Check number of rows of d (one fit per row) and
% whether z is a vector or matrix:

nrows = size(d, 1);


% Create the indices to call the rows of x. If it is a vector (1 row),
% indrx is a vector of 1's, otherwise it is a vector with its rows:
if isvector(x) || iscolumn(x)
    
    indrx = ones(1, nrows);
    
    if iscolumn(x)
        x = x';
    end
    
else
    
    if size(x, 1)~=nrows
        error('Arrays x and d are not consistent.')
    else
        indrx = 1:nrows;
    end
    
end


%% Check whether wnd and wndptslide were specified appropriately:

if ~exist('wnd', 'var')
    
    if exist('wndptslide', 'var')
        error('This does not make sense!')
    end
    
end


%% Create structure variable that goes
% into the least squares fit function:

% Constant:
imfWave.power = 0;

% Sinusoidal with frequency 1/prd:
imfWave.sine = 1/prd;


%% Do the fit:

% If there is no window, do simple least square fit for each row:
if ~exist('wnd', 'var')
    
    % Pre-allocate space for the wave-averaged variables:
    davg = NaN(nrows, 1);
    
    % Loop through rows of d:
    for i = 1:nrows

        % Fit:
        [~, wvavg_aux] = myleastsqrs(x(indrx(i), :), d(i, :), imfWave);

        % Assign the mean of the sinusoidal to the output variable:
        davg(i, 1) = wvavg_aux(1);
    end

    % Create the mean time output: 
    xavg = mean([x(1) x(end)]);
    
% If there is a window, use sliding harmonic fit to wave-average:
else
    
    [davg, xavg] = sliding_harmonicfit(x, d, wnd, wndptslide, ...
                                                    imfWave, [true false]);
    
end


                                                




