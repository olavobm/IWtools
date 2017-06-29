function [cp, cg, cn] = cteN_cpcg(freq, N, lat, D, nmd)
% [cp, cg, cn] = CTEN_CPCG(freq, N, f, D, nmd)
%
%   inputs:
%       - freq: wave frequency.
%       - N: buoyancy frequency.
%       - f: Coriolis parameter.
%       - D: fluid depth.
%       - nmd: mode number.
%
%   outputs:
%       - cp: phase speed (in m/s) of the mode n.
%       - cg: group velocity.
%       - cn: eigenspeed
%
% Olavo Badaro Marques, 28/Jun/2017.


%% Coriolis parameter:

f = gsw_f(lat);


%% Eigenspeed, phase speed and groud velocity of mode n:

cn = cteN_cn(freq, N, f, D, nmd);

[cp, cg] = cn2cpcg(cn, freq, lat);