function wvAng = iwChar(wvf, N, f0)
% wvAng = IWCHAR(wvf, N, f0)
%
%   inputs:
%       - wvf: wave frequency.
%       - N: buoyancy frequency.
%       - f0: Coriolis parameter
%
%   outputs:
%       - wvAng: first-quadrant angle (in radians)
%                of the internal wave characteristc.
%
% Returns the first-quadrant angle between the internal-wave
% group speed and the horizontal plane. This is given by the
% dispertion relationship of internal waves.
%
% Inputs can either be scalars or arrays. If more than one
% are (non-scalar) arrays, then they must have the same size.
%
% All inputs must have the same units and f0 <= wvf <= N. If this
% condition is not met, the output is NaN (in the appropriate
% locations, if inputs are arrays).
%
% Olavo Badaro Marques, 04/Mar/2017.


%% The dispersion relationship gives the tangent.

wvTan = abs( sqrt((wvf.^2 - f0.^2)./(N.^2 - wvf.^2)) );

wvAng = atan(wvTan);


%% See (if) where magnitude of the inputs
% is inconsistent with freely propagating
% internal waves:

lnonpropwaves = (f0 > wvf) | (wvf > N);

wvAng(lnonpropwaves) = NaN;