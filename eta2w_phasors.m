function wphasor = eta2w_phasors(wvfreq, etaphasor)
% wphasor = ETA2W_PHASORS(wvfreq, etaphasor)
%
%   inputs
%       - wvfreq: wave frequency in radians per second.
%       - etaphasor:
%
%   outputs
%       - wphasor:
%
% .... Linear approximation .... wave phasor with
% the convention of exp(-i*wvfreq*t)
%
% Olavo Badaro Marques, 24/Apr/2018.

wphasor = (-1i) .* wvfreq .* etaphasor;