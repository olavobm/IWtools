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
%       - slopecritc: criticality of the bottom.
%
% The criticality is the ratio between the (linear) bottom slope and the
% internal wave characteristic. The ranges of this parameters are:
%       * slopecritc >> 1: supercritical.
%       * slopecritc << 1: subcritical.
%       * slopecritc ~ 1: critical.
%
% A supercritical (subcritical) slope reflects
% the wave backwards (forward).
%
% See also: iwCharacteristic.m
%
% Olavo Badaro Marques, 24/Nov/2017.


%% Compute internal characteristic slope

iwRay = iwCharacteristic(wvfreq, N2, f0);


%% Ratio between bottom slope and wave characteristic

slopecritc = abs(bottomSlope) ./ iwRay;




