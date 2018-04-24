function [w, midtime] = eta2w(time, eta)
% [w, midtime] = ETA2W(time, eta)
%
%   inputs
%       - time: time
%       - eta: displacement
%
%   outputs
%       - w: d(eta)/dt
%       - midtime:
%
%
%
%
%
%
% Olavo Badaro Marques, 24/Apr/2018.


%%

if isvector(time) && ~isvector(eta)
	time = repmat(time, size(eta, 1), 1);
end

%
midtime = (time(:, 1:end-1) - time(:, 2:end)) ./ 2;


%%
w = (eta(:, 2:end)  - eta(:, 1:end-1)) ./ ...
    (time(:, 2:end) - time(:, 1:end-1));
