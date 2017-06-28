function wvhlen = cteN_wvhlen(freq, N, f, D, nmd)
% wvhlen = CTEN_WVHLEN(freq, N, f, D, nmd)
%
%   inputs:
%       - freq: wave frequency (in radians per second).
%       - N: buoyancy frequency.
%       - f: Coriolis parameter.
%       - D: fluid depth.
%       - nmd: mode number.
%
%   outputs:
%       - wvhlen: horizontal wavelength (in meters).
%
%
%
% Olavo Badaro Marques, 28/Jun/2017.


%% Quantized vertical wavenumber:

m_wvnmbr = nmd * pi / D;


%% Ratio of vertical with the horizontal wavenumber:

ratioWvNmbrs = sqrt((N.^2 - freq.^2) ./ (freq.^2 - f.^2));


%% Horizontal wavenumber magnitude:

h_wvnmbr = m_wvnmbr ./ ratioWvNmbrs;


%% Horizontal wavelength:

wvhlen = 2*pi ./ h_wvnmbr;


