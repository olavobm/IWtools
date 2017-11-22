function vurto = ellipseLat0(lat0, wvfreq)
% vurto = ELLIPSELAT0(lat0, wvfreq)
%
%   inputs
%       - lat0: latitude (in degrees).
%       - wvfreq: wave frequency (in cycles per day).
%
%   outputs
%       - vurto: ratio magnitude between amplitude of the
%                velocity components perpendicular and parallel
%                to wave propoagation.
%
%
% Dependecies: gsw_f.m.
%
% Olavo Badaro Marques, 21/Nov/2017.


%%

f0 = gsw_f(lat0);


%%

if isa(wvfreq, 'char')
    wvfreq = tidalFreq(wvfreq);
end

wvfreq = (2*pi/(24*3600))* wvfreq;



%%

vurto = abs(f0) ./ wvfreq;


