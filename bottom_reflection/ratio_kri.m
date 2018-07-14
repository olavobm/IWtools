function krirto = ratio_kri(wvchar_angle, bottom_angle)
% krirto = RATIO_KRI(wvchar_angle, bottom_angle)
%
%   inputs
%       - wvchar_angle: angle (in radians) of the wave
%                       characteristic with the horizontal.
%       - bottom_angle: angle (in radians) of the linearly
%                       sloping bottom.
%
%   outputs
%       - krirto: ratio between reflected and incident horizontal
%                 wavenumbers of the waves on a linearly sloping
%                 bottom.
%
% The equation for krirto comes from the condition of no normal
% flow through the boundary (first derived by Phillips, 1977).
%
% See also: iwChar.m
%
% Olavo Badaro Marques, 14/Jul/2018.


%%

%
add_angles = wvchar_angle + bottom_angle;

%
absdiff_angles = abs(wvchar_angle - bottom_angle);

%
krirto = sin(add_angles) ./ sin(absdiff_angles);

