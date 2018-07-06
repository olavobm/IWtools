function [pp, psurf] = eta2pp(z, eta, N2, zN2, rho0)
%% [pp, psurf] = ETA2PP(z, eta, N2, zN2, rho0)
%
%  inputs
%    - z: vector of the data points depths in meters. Depth is greater
%         than 0 and should be specified in ascending order.
%    - eta: vector or matrix of the isopycnal displacement.
%    - N2: background buoyancy frequency squared. If it is a
%          vector, it is applied everywhen in eta.
%    - zN2 (optional): depth where N2 is specified. If not given, code
%                      assumes N2 is given at the same depth as eta.
%    - rho0 (optional): reference potential density (default is 1025).
%
%  outputs
%    - pp: pressure perturbartion.
%    - psurf: baroclinic surface pressure.
%
% ETA2PP computes the pressure perturbation profile due to hydrostatic
% internal wave (hydrostatic valid when wavefrequency^2 << N^2).
% The calculation assumes the depth integral of pp is zero.
%
% The calculation is equation (1) from Kunze et al. (2002): Internal
% Waves in Monterey Submarine Canyon.
%
% TO DO:
%   - Fix size/dimensions of inputs/outputs
%   - Use different boundary conditions???
%   - Clean the code.
%
% Olavo Badaro Marques, 27/Oct/2016.


%% Define default value for reference density:

if ~exist('rho0', 'var')
    rho0 = 1025;
end


%% If input zN2 is given, then interpolate N2 to the same depths as eta:

if ~exist('zN2', 'var')
    
else
    
    N2 = interp1(zN2, N2, z);
    
end


%% Maybe resize an eta cube array to a 2D matrix in
% order to simplify the code...


%% Check dimensions of eta and N2:

% If eta is only a vector, make
% sure it is a column vector:
if isrow(eta)
	eta = eta'; 
end

% Arrange N2 array:
if isvector(N2)
    
    % Make sure N2 is a column vector:
	N2 = N2(:);
    
    % Make it a matrix with the same size as the displacement:
    N2 = repmat(N2, 1, size(eta, 2));
    
end


%% Compute the perturbation density multiplied
% by the gravity acceleration:

rhopg = rho0 .* eta .* N2;


%% 

% Pressure perturbation (IMPROVE FOR UPPER NaNs!!!):
pptop = zeros(1, size(rhopg, 2));

pp = NaN(size(rhopg));
psurf = NaN(1, size(rhopg, 2));

%
ncols = size(rhopg, 2);
allcols = 1:ncols;
allcols = allcols';  % make it column vector

% Get column subscripts of where there are NaNs:
[~, cnan] = ind2sub(size(rhopg), find(isnan(rhopg)));
cNanNoRep = unique(cnan);

cGoodData = setdiff(allcols, cNanNoRep);


%% Cumulatively integrate the density perturbation:

intRho = NaN(size(rhopg));

% Deal with all profiles that do not have NaNs
if ~isempty(cGoodData)
    intRho(:, cGoodData) = cumtrapz(z, rhopg(:, cGoodData));
    
    Dzrange = z(end)-z(1);
    
    % Boundary condition:
    psurf(cGoodData) = pptop(cGoodData) - (1/Dzrange) .* trapz(z, intRho(:, cGoodData));
    
    % Turn psurf into a matrix to add to intRho:
	psurf = repmat(psurf, size(rhopg, 1), 1);
    
    
    % Compute the pressure perturbation
    pp(:, cGoodData) = psurf(:, cGoodData) + intRho(:, cGoodData);
    
end

% For the profiles that have NaNs, do one by one
if ~isempty(cNanNoRep)
    
    Dzrange = NaN(1, size(rhopg(:, cNanNoRep), 2));
    
    for i = 1:size(rhopg(:, cNanNoRep), 2)
        
        %
        lgood = ~isnan(rhopg(:, cNanNoRep(i)));
        zgood = z(lgood);
        
        %
        if any(lgood)
            intRho(lgood, cNanNoRep(i)) = cumtrapz(zgood, rhopg(lgood, cNanNoRep(i)));
            try
            Dzrange(i) = zgood(end)-zgood(1);
            catch
                keyboard
            end
            % Boundary condition:
            psurf(1, cNanNoRep(i)) = pptop(cNanNoRep(i)) - (1/Dzrange(i)) .* trapz(zgood, intRho(lgood, cNanNoRep(i)));

            % NEED TO DEAL WITH NANS PROPERLY!!!!

            % Compute the pressure perturbation (no need to use repmat
            % since we are dealing with profiles individually):
            pp(lgood, cNanNoRep(i)) = psurf(1, cNanNoRep(i)) + intRho(lgood, cNanNoRep(i));
        end
    end
    
end



