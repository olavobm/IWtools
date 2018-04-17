function y = phasor2timeseries(t, freq, x)
% y = PHASOR2TIMESERIES(t, freq, x)
%
%   inputs
%       - t: 
%       - freq: 
%       - x: the phasor (a complex number).
%
%   outputs
%       - y:
%
% PHASOR2TIMESERIES.m computes a timeseries from the phasor
% (i.e. complex amplitude) x. The convention in this code
% is that y = Re{x*exp[i*freq*t]}.
%
% y = abs(x) * cos(freq*t + .....) % {atan2(imag(x), real(x))]}. Correct this.
% 
%
% Note that wave solutions are often expressed in terms
% of exp(-freq*t), so that in this case the input freq
% should be negative.
%
% Olavo Badaro Marques, 17/Apr/2018.


y = real(x).*cos(freq.*t) - imag(x).*sin(freq.*t);

