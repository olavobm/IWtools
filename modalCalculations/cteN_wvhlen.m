function wvhlen = cteN_wvhlen(wvfreq, N, lat0, D, nmd)
% wvhlen = CTEN_WVHLEN(wvfreq, N, f, D, nmd)
%
%   inputs
%       - wvfreq: wave frequency (in cycles per day).
%       - N: buoyancy frequency.
%       - lat0: latitude (in degrees).
%       - D: water depth.
%       - nmd: mode number.
%
%   outputs
%       - wvhlen: horizontal wavelength (in meters).
%
%
%
% Olavo Badaro Marques, 28/Jun/2017.


%% Convert from cycles per day to radians per second

wvfreq = wvfreq * 2*pi / (24*3600);


%% Compute Coriolis parameter

f0 = gsw_f(lat0);


%% Quantized vertical wavenumber:

m_wvnmbr = nmd * pi / D;


%% Ratio of vertical with the horizontal wavenumber:

ratioWvNmbrs = sqrt((N.^2 - wvfreq.^2) ./ (wvfreq.^2 - f0.^2));


%% Horizontal wavenumber magnitude:

h_wvnmbr = m_wvnmbr ./ ratioWvNmbrs;


%% Horizontal wavelength:

wvhlen = 2*pi ./ h_wvnmbr;


