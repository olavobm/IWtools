function oceanModes = oceanVmodes(z, H, N2, nmds)
%% oceanModes = OCEANVMODES(z, H, N2, nmds)
% 
%   inputs
%       - z: depth grid vector where N2 is specified on (can be irregular).
%            The first (last) element MUST BE greater (less) than 0 (H),
%            to emphasize that N2 at z = 0 (z = H) is irrelevant for the
%            calculation.
%       - H: total water depth.
%       - N2: vector with buoyancy frequency squared. NO NaNs allowed.
%       - nmds: number of modes to extract (from 0 to nmds).
%
%   outputs
%       - oceanModes: structure variable with the following field
%           * z: depth vector associated with the modes.
%           * N2: buoyancy frequency squared (same as input).
%           * nmode: number of the mode (columns of the eigenfunctions).
%           * Heigfcn: matrix of size (length(z)+2, nmds), where
%                      each column is an eigenfunction of pressure
%                      (and horizontal velocity). Increasing mode with
%                      increasing column index.
%           * Veigfcn: same as Heigfcn, but for the eigenfunctions of
%                      buoyancy (and vertical velocity).
%           * ce: eigenspeeds.
%
% TO DO:
%       output - He: equivalent depth... I need rotation to correct g.
%
%       input - lmode0: Set true (default) to include mode-0
%                       set false to exclude mode-0
%
% The eigenmode outputs are normalized such that each mode is bounded
% by -1 and +1, which makes them non-orthogonal (IN GENERAL?? WHAT ABOUT
% CONSTANT STRATIFICATION???). Note that their number
% of rows is equal to length(z)+2, because the input z does not include
% the surface and bottom whereas the output modes include values there.
%
% This function computes the inertia-gravity waves vertical standing
% modes for an "arbitrary" (??????) buoyancy frequency squared, N2. It
% is assumed that the frequency of the wave is much smaller than the
% buoyancy frequency, which implIes a hydrostatic approximation. The
% vertical orthogonal modes associated with the vertical velocity are
% solutions of the equation:
%
%   V_zz + ev*N2*V = 0,
%
% with boundary conditions (use real or rigid-lid????
%                           should have both possibilities????):
%
%   V(H) = 0    and    Vz(0) - g*ev*V(0) = 0.
%
% Physically, the eigenvalues ev are the inverse of the eigenspeeds
% squared. The eigenspeed is the square root of the equivalent depth
% times the acceleration of gravity. It is also the geometric mean
% of the phase speed and group velocity of each mode. If the planetary
% rotation is zero, all these speeds are numerically equal.
% 
% Medium matrix (A). This is the one that multiplies V in the
% equation and contains information about the vertical variation
% of the medium (i.e. stratification) where the wave is propagating.
%
% - IS IT ONLY SYMMETRICAL WHEN DZ(1) IS THE SAME AS THE SURFACE GAP??????
% - NON-NORMALIZED EIGENFUNCTIONS?
% - WHEN SOLVING THE EIGENVALUE PROBLEM, CAN I SOLVE FOR 0TH MODE??????
%
% Note sw_vmodes computes the eigenfunctions in a different, clever way.
%
% CAVEAT: Solves an eigenvalue equation, which means it will create an mxm
% matrix and solve it.  Don't put in a 4000m drop with data every meter;
% smooth and decimate first!...... UPDATE THIS CAVEAT COMMENT
%
% Functions called by OCEANVMODES: EignProblemD2 and centeredDeriv.
%
% based on MODES.m by Sam Kelly.
% Olavo Badaro Marques, 31/10/2016.


%% Assign N2 to another variable, which is used in the calculation.
% In the case where there are surface/bottom values, which do not go
% into the calculation, this assignment allows the original N2 to
% be assigned to the output variable

N2diag = N2;


%% Check inputs

if ~isvector(z)
    error('Input z is not a vector!')
else
    if isrow(z)
        z = z(:);
    end
    
    if z(1)<0
        error(['The first element of the depth vector must ' ...
               'be greater or equal than 0 (the surface).'])
    end
    
    if z(end)>H
        error(['The last element of the depth vector must '   ...
               'be less or equal than the water depth.'])
    end
    
    % If first (last) depth levels are top (bottom), remove
    % them. Afterall, they are not used in the calculation
    if z(1)==0;    z = z(2:end);   N2diag = N2diag(2:end);    end
    if z(end)==H;  z = z(1:end-1); N2diag = N2diag(1:end-1);  end
    
end

nz = length(z);


%% If not specified as input, choose number of modes. Since the
% maximum number of modes that can be computed is limited by
% the number of observations, set a fixed max number of modes
% in case there are too many data point

% FORMALLY, WHAT IS THE CONSTRAIN IN THE NUMBER OF MODES???????????
% for n = 3, it's fine to get 3 modes. For n = 35, it only works up
% to 33 modes.

maxnmds = 100;

if ~exist('nmds', 'var')
    if nz > maxnmds
        nmds = maxnmds;
    else
        nmds = nz - 2; % I think it can be n
    end
else
    % make sure input nmds is no more than we can compute.
end


%% Check whether z is regularly spaced (in which case
% the derivative matrix is symmetric and the solution
% of the eigenvalue problem can be optimized)


% EVEN IF IT IS REGULAR, IT IS PROBABLY NEVER GOING TO BE REGULAR
% TO THE BOTTOM DEPTH!!!!!!! DEAL WITH THAT!!!!

dz = z(2:end) - z(1:end-1);

% If z is regular -> symmetric matrix:
lsym = all(dz == dz(1));   % true if z is regular

% if all(dz == dz(1))    
%     
% end


%% Second derivative matrix for a general position vector,
% assuming homogeneous Dirichlet boundary conditions

% Add surface and bottom to vector z, where we have boundary conditions:
zsb = [0; z; H];
n = nz + 2;

% Call another function to make the derivative matrix:
D2 = EignProblemD2(zsb);


%% Medium matrix A

A = diag(-N2diag);


%% Solve generalized eigenvalue problem

opts.issym = lsym;
opts.isreal = 1;

[Veigfcn, eigvals2] = eigs(D2, A, nmds, 'SM', opts);


%% Sort modes by eigenspeed

% Get the diagonal values with the eigenvalues:
eigvals2 = diag(eigvals2);
eigvals2(eigvals2<0) = Inf;  % when would this happen??? I think this is
                             % that evanescent modes have zero eigenspeed

% Compute the eigenspeed:
ce = 1 ./ sqrt(eigvals2);

% Sort eigenspeeds and eigenfunctions from highest to lowest
% (actually, the option 'SM' in eigs alreay returns sorted output):
[ce, ind] = sort(ce, 'descend');
Veigfcn = Veigfcn(:, ind);

% Add the boundary values at the surface (rigid lid) and
% bottom, which satisfy the boundary conditions:
Veigfcn = [zeros([1, nmds]); ...
               Veigfcn    ; ...
           zeros([1, nmds])];

       
%% Scale Veigfcn appropriately
%
% ----- SHOULD FIX THIS FOR IRREGULAR Z -----

ffactor = sum(repmat([0; N2diag; 0], 1, nmds) .* Veigfcn.^2 .* mean(dz) , 1);
ffactor = sqrt(ffactor);

Veigfcn = Veigfcn ./ repmat(ffactor, size(Veigfcn, 1), 1);

% % % This integral for each mode is 1
% % sum(repmat([0; N2diag; 0], 1, nmds) .* bla.^2 .* mean(dz) , 1)


%% Compute velocity/pressure modes by taking the derivative
% of displacement modes. Then scale them appropriately.
%
% As far as I understand, it should not be necessary to scale
% Heigfcn. My guess is that Heigfcn is "approximately" scaled
% (integral close to 1), but maybe it is something else.
       
% Take derivative to get U and P structure (forward/backward
% differences at the edges and centered differences in the interior):
Heigfcn = centeredDeriv(zsb, Veigfcn);

% Scale Heigfcn
ffactor = sum(Heigfcn.^2 .* mean(dz) , 1);
ffactor = sqrt(ffactor);

Heigfcn = Heigfcn ./ repmat(ffactor, size(Heigfcn, 1), 1);

% % % This integral for each mode is 1
% % sum(Heigfcn.^2 .* mean(dz) , 1)


% ------------------------------------------------------------------------
% % %% Normalization of the eigenfunctions
% % % such that the are bounded by -1 and +1
% % 
% % Heigmax = max(abs(Heigfcn));
% % Veigmax = max(abs(Veigfcn));
% % 
% % Heigfcn = Heigfcn ./ repmat(Heigmax, [n 1]);
% % Veigfcn = Veigfcn ./ repmat(Veigmax, [n 1]);
% ------------------------------------------------------------------------


%% Multiply by -1 the modes that have negative
% surface Heigfcn (simply for aesthetic puporses)

lcols = Heigfcn(1, :) < 0;
Heigfcn(:, lcols) = -Heigfcn(:, lcols);
Veigfcn(:, lcols) = -Veigfcn(:, lcols);


%% Add 0th mode

Heigfcn = [ones(n, 1), Heigfcn];
Veigfcn = [(H-zsb)./H, Veigfcn];
    
ce = [sqrt(9.81*H); ce];


%% Assign output variable

oceanModes.z = zsb;
oceanModes.N2 = N2;
oceanModes.nmode = 0 : nmds;
oceanModes.Hfcn = Heigfcn;
oceanModes.Vfcn = Veigfcn;
oceanModes.ce = ce;

