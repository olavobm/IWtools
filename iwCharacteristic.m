function raySlope = iwCharacteristic(wvfreq, N2, f0)
% raySlope = IWCHARACTERISTIC(wvfreq, N2, f0)
%
%   inputs
%       - wvfreq: wave frequency (in cycles per day).
%       - N2: buoyancy frequency squared (in radians per second squared).
%       - f0: Coriolis parameter (in radians per second).
%
%   outputs
%       - raySlope: tangent of the angle the wave characteristic
%                   makes with the horizontal.
%
% Olavo Badaro Marques, 24/Nov/2017.


%% Convert wave frequency from cycles per day to radians per second

wvfreq = (2*pi/(24*3600))* wvfreq;


%% Check if  N2 >= wvfreq^2 >= f0^2 (condition for a free wave)

if ~( (sqrt(N2) >= wvfreq) && (wvfreq >= f0) )
    error('Wave frequency is not bounded by N2 and f0.')
end


%% Compute (positive) tangent of the angle between the
% internal wave characteristic and the horizontal

raySlope = abs( sqrt((wvfreq.^2 - f0.^2) ./ (N2 - wvfreq.^2)) );