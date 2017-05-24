function [pp, pbc] = presIW_fromEta(z, eta, N2, rho0)
%% pp = PRESIW_FROMETA(z, eta)
%
%  inputs:
%    - z: vector of the data points depths in meters. Depth is greater
%         than 0 and should be specified in ascending order.
%    - eta: vector or matrix of the isopycnal displacement.
%    - N2: background buoyancy frequency squared. If it is a
%          vector, it is applied everywhen in eta.
%    - rho0 (optional): reference potential density (default is 1025).
%
%  outputs:
%    - pp: pressure perturbartion.git 
%    - pbc:
%
% The function PRESIW_FROMETA computes the pressure perturbation profile
% due to hydrostatic internal wave (hydrostatic in the sense that
% wavefrequency^2 << N^2).
%
% --------------------------   NaNs????  --------------------------
%
% The calculation assumes the depth integral of pp is zero. ??????????
%
% Olavo Badaro Marques, 27/Oct/2016.


%%

if ~exist('rho0', 'var')
    rho0 = 1025;
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

%% Create the boundary condition according to normal modes (N2)


%% Compute the perturbation density multiplied
% by the gravity acceleration:

rhopg = rho0 .* eta .* N2;


%% .........

% Cumulatively integrate the perturbation density:
% and the boundary condition:

% Pressure perturbation (IMPROVE FOR UPPER NaNs!!!):
pptop = zeros(1, size(rhopg, 2));

pp = NaN(size(rhopg));

% also prea-allocate space for pbc

% -------------------------------------------------------------------------
% ACTUALLY, IT IS BETTER TO DO THE PROFILE OPTIMIZATION THING, AFTERALL
% IF THERE ARE NO NANS, ALL PROFILES WILL GO ON THE SET WHERE WE CAN DO
% RIGHT AWAY AND THE OTHER SET WILL BE EMPTY.


ncols = size(rhopg, 2);
allcols = 1:ncols;
allcols = allcols';  % make it column vector

% Get column subscripts of where there are NaNs:
[~, cnan] = ind2sub(size(rhopg), find(isnan(rhopg)));
cNanNoRep = unique(cnan);

cGoodData = setdiff(allcols, cNanNoRep);

% -------
intRho = NaN(size(rhopg));

if ~isempty(cGoodData)
    intRho(:, cGoodData) = cumtrapz(z, rhopg(:, cGoodData));
    
    Dzrange = z(end)-z(1);
    
    % Boundary condition:
    pbc = pptop(cGoodData) - (1/Dzrange) .* trapz(z, intRho(:, cGoodData));
    
    % Make a pbc matrix (here, it is almost always the
    % case that pbc is a matrix and not a vector):
	pbc = repmat(pbc, size(rhopg, 1), 1);
    
    % Compute the pressure perturbation
    pp(:, cGoodData) = pbc + intRho(:, cGoodData);
end

if ~isempty(cNanNoRep)
    
    Dzrange = NaN(1, size(rhopg(:, cNanNoRep), 2));
    
    for i = 1:size(rhopg(:, cNanNoRep), 2)
        
        lgood = ~isnan(rhopg(:, cNanNoRep(i)));
        zgood = z(lgood);
        intRho(lgood, cNanNoRep(i)) = cumtrapz(zgood, rhopg(lgood, cNanNoRep(i)));
        
        Dzrange(i) = zgood(end)-zgood(1);
        
        % Boundary condition:
        pbc = pptop(cNanNoRep(i)) - (1/Dzrange(i)) .* trapz(zgood, intRho(lgood, cNanNoRep(i)));

        % NEED TO DEAL WITH NANS PROPERLY!!!!
        
        % Compute the pressure perturbation (no need to use repmat
        % since we are dealing with profiles individually):
%         try
        pp(lgood, cNanNoRep(i)) = pbc + intRho(lgood, cNanNoRep(i));
%         catch
%             keyboard
%         end
        
    end
    
end




% % % % %% In case there are NaNs, we should deal with each profile individually
% % % % % (or at least all profiles with NaNs inidividually):
% % % % if any(isnan(rhopg(:)))
% % % %     
% % % %     
% % % %     
% % % % %% If there are NO NaNs, we can deal with all profiles at once:   
% % % % else
% % % %             
% % % %     intRho = cumtrapz(z, rhopg);
% % % %     
% % % %     Dzrange = z(end)-z(1);
% % % %     
% % % %     % Boundary condition:
% % % %     pbc = pptop - (1/Dzrange) .* trapz(z, intRho);
% % % %     
% % % %     % Make a pbc matrix (here, it is almost always the
% % % %     % case that pbc is a matrix and not a vector):
% % % % 	pbc = repmat(pbc, size(rhopg, 1), 1);
% % % %     
% % % %     % Compute the pressure perturbation
% % % %     pp = pbc + intRho;
% % % %     
% % % % end


