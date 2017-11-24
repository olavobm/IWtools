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

raySlope = abs( sqrt((wvfreq^2 - f0^2)/(N2 - wvfreq^2)) );