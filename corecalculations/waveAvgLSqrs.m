function [davg, xavg] = waveAvgLSqrs(x, d, prd, porder, wnd, wndptslide)
% [davg, xavg] = WAVEAVGLSQRS(x, d, prd, wnd, wndptslide)
% 
%   inputs:
%       - x: vector or matrix.
%       - d: variable to be wave-averaged (each of its rows).
%       - prd: period of the sinusoidal (or vector of periods), in
%              the same units as x.
%       - porder (optional): order of polynomials to fit, rather than just
%                            fitting a constant, such that davg is the mean
%                            of the polynomial fit.
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

%
xncols = size(x, 2);


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

% Sinusoidal(s) with frequency(ies) 1./prd:
imfWave.sine = 1./prd;

%
if ~exist('porder', 'var') || isempty(porder)
else
    if iscolumn(porder)
        porder = porder';
    end
    imfWave.power = [imfWave.power porder];
end


%% Do the fit:

% If there is no window, do simple least square fit for each row:
if ~exist('wnd', 'var')
    
    % Pre-allocate space for the wave-averaged variables:
    davg = NaN(nrows, 1);
    
    % Loop through rows of d:
    for i = 1:nrows

        % Fit:
        [~, wvavg_aux] = myleastsqrs(x(indrx(i), :), d(i, :), imfWave);

        %
        if ~exist('porder', 'var')        
            % Just assign the constant from the fit to the output variable:
            davg(i, 1) = wvavg_aux(1);
        else
            % In this case, reconstruct the polynomial part of the
            % fit on a higher resolution grid and take the mean:
            nauxgrid = 2*xncols;
            x2G = linspace(x(indrx(i), 1), x(indrx(i), end), nauxgrid);
            x2G = x2G';
            x2G = repmat(x2G, 1, length(imfWave.power));
            
            Gaux = x2G .^ repmat(imfWave.power, nauxgrid, 1);
            
            davg(i, 1) = mean(Gaux*wvavg_aux(1:length(imfWave.power)));
            
            % if x is just a vector, the creation of Gaux can be done
            % before the loop.
        end
        
    end

    % Create the mean time output: 
    xavg = mean([x(1) x(end)]);
    
% If there is a window, use sliding harmonic fit to wave-average:
else

    [davg, xavg] = sliding_harmonicfit(x, d, wnd, wndptslide, ...
                                            imfWave, [true, false, false]);
    
end


                                                




