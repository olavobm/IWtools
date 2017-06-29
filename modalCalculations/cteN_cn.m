function cn = cteN_cn(freq, N, f, D, nmd)
% cn = CTEN_CN(freq, N, f, D, nmd)
%
%   inputs:
%       - freq: wave frequency.
%       - N: buoyancy frequency.
%       - f: Coriolis parameter.
%       - D: fluid depth.
%       - nmd: mode number.
%
%   outputs:
%       - cn: eigenspeed of the mode n.
%
% All frequencies must be given in radians per second.
%
% Olavo Badaro Marques, 28/Jun/2017.

iwRatio = (freq.^2 - f.^2) ./ (N.^2 - freq.^2);

cn_num = N.^2 - f.^2;
cn_den = ((nmd .* pi / D).^2) .* ( iwRatio + 1);

cn = cn_num ./ cn_den;