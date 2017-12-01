function vurto = ellipseLat0(lat0, wvfreq)
% vurto = ELLIPSELAT0(lat0, wvfreq)
%
%   inputs
%       - lat0: latitude (in degrees).
%       - wvfreq: wave frequency (in cycles per day).
%
%   outputs
%       - vurto: magnitude of the ratio between the amplitudes of
%                the velocity components perpendicular and parallel
%                to wave propoagation.
%
% From the equations of motion, inertia-gravity waves (barotropic or
% each baroclinic mode) have v/u = -i*f0/wvfreq -- which describes
% clockwise particle trajectories in the northern hemisphere. In this
% notation, u is horizontal velocity component in the direction of
% wave propagation.
%
% Dependecies: gsw_f.m.
%
% Olavo Badaro Marques, 21/Nov/2017.


%% Computes the Coriolis parameter

f0 = gsw_f(lat0);


%% Convert wave frequency from cycles per day to radians per second

wvfreq = (2*pi/(24*3600))* wvfreq;


%% Compute the magnitude of the ratio between
% amplitude of the velocity components

vurto = abs(f0) ./ wvfreq;


