function slopecritc = slopeCriticality(bottomSlope, wvfreq, N2, f0)
% slopecritc = SLOPECRITICALITY(bottomSlope, wvfreq, N2, f0)
%
%   inputs
%       - bottomSlope: tangent of the angle between the
%                      bottom and the horizontal.
%       - wvfreq: wave frequency (in cycles per day).
%       - N2: buoyancy frequency squared (in radians per second squared).
%       - f0: Coriolis parameter (in radians per second).
%
%   outputs
%       - slopecritc: ratio between bottom slope and
%                     internal wave characteristic.
%
% The criticality is the ratio .... The criticality is classified as
%       * slopecritc >> 1: supercritical.
%       * slopecritc << 1: subcritical.
%       * slopecritc ~ 1: critical.
%
%
% Dependecies: iwCharacteristic.m
%
% Olavo Badaro Marques, 24/Nov/2017.


%% Compute internal characteristic slope

iwRay = iwCharacteristic(wvfreq, N2, f0);


%% Ratio between bottom slope and wave characteristic

slopecritc = abs(bottomSlope) ./ iwRay;




