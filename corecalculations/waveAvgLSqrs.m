function [davg, xavg] = waveAvgLSqrs(x, d, prd) %, wnd)
% [davg, xavg] = WAVEAVGLSQRS(x, d, prd, wndw)
% 
%   inputs:
%       - x: time vector correspondent to each column of d
%       - d: data of the observations
%       - prd: period of the sinusoidal 
%       - wnd [optional]: window to 
%
%   outputs:
%       - davg:
%       - xavg:
%   MAYBE OUTPUT THE FIT AS WELL
%
% Use least squares to compute the wave-average (davg) of d. The
% coefficients of a sinusoidal of period prd plus a constant are
% estimated from least squares and the constant is the wave-averaged
% output is the constant.
% The fit is computed by the function myleastsqrs.m.
%
% - d can be either vector or matrix.
% - units of x and prd must be consistent among each other.
% - wnd is a window. This can also be used as filtering as well.
%
% If no window is specified, the output is just a number or vector with
% the wave-averaged quantities.
%
% Olavo Badaro Marques, 17/08/2016.

% MAYBE I SHOULD CALL THE FUNTION SLIDING_HARMONCFIT INSIDE THIS FUNCTION

%% Check number of rows of d. We will do one fit per row:
%
% COLUMN VECTOR???

Nfits = size(d, 1);


%% Check window input:


%% Create structure variable that goes
% into the least squares fit function:

% Constant:
imfWave.power = 0;

% Sinusoidal with frequency 1/prd:
imfWave.sine = 1/prd;


%% Pre-allocate space for wave-averaged quantity:

davg = NaN(Nfits, 1); % IT COULD ALSO BE A MATRIX IF WE SPECIFY A WINDOW!!!
% xavg = NaN(1, 1);     % "    "     "  BE A VECTOR IF WE    "    "    "!!!


%% Do the fit -- outside for loop goes through the windows
% of the data and inside for loop through the rows of it:

% Loop through rows of d:
for i = 1:Nfits
    
    % Fit:
    [~, wvavg_aux] = myleastsqrs(x, d(i, :), imfWave);
    
    % Assign the mean of the sinusoidal to the output variable:
    davg(i, 1) = wvavg_aux(1);
    
end

% Create the mean time output: 
xavg = mean(x);

