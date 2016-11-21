function [MdsProf, MdsAmp, VarMds] = modalDecomp(nmds, zN2, N2, D, varStruct, N2tandwnd)
% MODALDECOMPOSITION(zN2, N2, D, zH, varH, zV, var, N2window)
%
% varStruct????
%
%   inputs:
%       - nmds: number of BAROCLINIC(???) modes to fit (FIX THIS PURPOSE!!!)
%       - zN2:
%       - N2: vector or matrix (be aware that this would ideally not
%                               contain NaNs....).
%       - D: water depth (positive number).
%       - varStruct: structure variable with (possible???) fields:
%               * z:
%               * varH
%               * varV
%               * (optional and fancy)
%               *   zH: cell array correspondent to varH.
%               *   zV: cell array correspondent to varV.
%
%       - N2tandwnd (optional): cell array, whose (1) first element is the
%                               time associated with the columns of N2 and
%                               the data and (2) the second is the window
%                               width, in the time units, for averaging N2.
%
%   outputs:
%       - MdsProf: struct variable with modal structures.
%       - MdsAmp: struct variable with modal amplitudes.
%       - VarMds: struct variable with projection of the variables
%                 onto modes.
%
% Function MODALDECOMPOSITION does ....
%
% HOW DO I DEAL WITH THE 0TH MODE OF VERTICAL-TYPE EIGENFUNCTIONS???
% I WOULD HAVE TO CHECK WHETHER THE VARIABLE HAS DEPTH-INTEGRAL ~= 0
%
% If N2 is a matrix and no window is specified, then normal modes are
% computed for every column of N2. The part that may take a long time
% 
% THE REALLY AWESOME THING, WOULD BE TO INCORPORATE MATRIX INPUT INTO
% MYLEASTSQRS FUNCTION, SUCH THAT MODAL PROJECTION CAN BE DONE RIGHT
% AWAY FOR TIME-CONSTANT N2!!! MAKE SURE IT WORKS..... which I can easily
% test (for xgd matrix) with:
%                   m = ( G' * G) \ ( G' * xgd);
%
% -------------------------------------------------------------------------
% Does it make faster/slower if I rearrange the parenthesis above???
% IN WHICH CASE, I MUST NOT HAVE NANS OR HAVE NANS AT THE SAME DEPTHS FOR
% ALL TIMES
% -------------------------------------------------------------------------
%
% Olavo Badaro Marques, 18/Nov/2016.


%% Preliminaries:

nN2profs = size(N2, 2);


%% Get the names of the variables:

if isfield(varStruct, 'varH')    
    namesVarH = fieldnames(varStruct.varH);
else
    namesVarH = {};
end

if isfield(varStruct, 'varV')
    namesVarV = fieldnames(varStruct.varV);
else
    namesVarV = {};
end

namesallvar = [namesVarH; namesVarV];
lforHvar = [true(length(namesVarH), 1); false(length(namesVarV), 1)];
typesallvar(lforHvar) = 'varH';   % string to indicate appropriate
typesallvar(~lforHvar) = 'varV';  % structure field in a dynamical "way"
                                  % (get correct terminology, ???)

% Number of profiles we have to fit:
nprofs = size(varStruct.(namesallvar(1)), 2);
% FOR NOW (AND PROBABLY FOREVER), LET'S TAKE ALL THE DATA GRIDDED IN TIME.


%% Check whether N2tandwnd was specified and determine
% how many chuncks of data will be used:

N2time = N2tandwnd{1};
N2wndw = N2tandwnd{2};

N2halfwndw = round(N2wndw/2);

% I SHOULD LOOK AT SPECTRA FUNCTION/SARAH'S CODE DEFINE
% THE OVERLAP AND HOW TO SLIDE THE WINDOW

% AND CREATE SOME SORT OF N2 MATRIX AND INDICES TO IDENTIFY THE WINDOW
% WHERE IT APPLIES!!!!

if size(N2, 2)==1
    % For the case of an individual profile:
    indN2col = ones(1, nprofs);
    ndiffMds = 1;
else
    indN2col = 1:nprofs;
    ndiffMds = nprofs;
end
% THE WINDOWED CASE IS MISSING FROM THE IF ABOVE!!!


% MdsProf, MdsAmp


%% Fit modes:

% Loop through all variable names and pre-allocate
% space for amplitude variables:
for i = 1:length(namesallvar)
    
    if lforHvar(i)
        nrows_aux = nmds + 1;
    else
        nrows_aux = nmds;
    end
    
    MdsAmp.(namesallvar(i)) = NaN(nrows_aux, nprofs);
    
end

% veczmds = [MMP.z; MMP.depth];

lN2toModes = true;

% I SHOULD LOOP THROUGH WINDOWS, GIVEN THERE ARE NO NANS AND I CAN SOLVE
% THE LEAST SQUARES EFFICIENTLY. I DON'T KNOW WHAT SHOULD I LOOP THROUGH!!! 
for i = 1:ndiffMds
    
    % If logical variable is true, (IF LOOPING THROUGH ndiffMds, if is
    % useless)
%     if lN2toModes
        [Hmodes, Vmodes, ~] = inertGravVmodes(zN2, D, ...
                                              N2(:, indN2col(i)), nmds);
%     end
    
    if ndiffMds==nprofs
        % BRUTE FORCE - FIT FOR EVERY PROFILE
    else
        
        
        
        
        
    end
    % I should subset 
    
	eta_mdsAmp(:, i) = fitVmodes(zlvs(2:end-1), MMP.etaM2(:, i), ...
                                 Vmodes(:, 2:end), veczmds);
                             
	u_mdsAmp(:, i) = fitVmodes(MMP.z, MMP.uM2(:, i), ...
                               Hmodes, veczmds);
                           
	v_mdsAmp(:, i) = fitVmodes(MMP.z, MMP.vM2(:, i), ...
                               Hmodes, veczmds);
    
    % If it takes too long, print progress to the screen:
    if mod(i, 200)==0
        disp(['Fitting modes to column ' num2str(i) ' of ' num2str(length(MMP.yday)) ''])
    end
    
      
end

Vmdsinterp = NaN(length(MMP.z), nmds);
Hmdsinterp = NaN(length(MMP.z), nmds+1);

% ------------------------------------------------------
% THE ASSIGNMENT BELOW IS WRONG FOR VMODES!!!
% I'M TAKING THE 0TH VERTICAL MODE FOR COMPUTING
% 1ST MODE ETA AND PRESSURE.
keyboard
% ------------------------------------------------------

for i = 1:(nmds+1)
    
    if i<(nmds+1)
        Vmdsinterp(:, i) = interp1(veczmds, Vmodes(:, i), MMP.z);
        eval(['etamds' num2str(i) ' = Vmdsinterp(:, i) * eta_mdsAmp(i, :);']); 
        eval(['pmds' num2str(i) ' = presIW_fromEta(MMP.z, etamds' num2str(i) ', MMP.n2extrap);']);
    end
  
    Hmdsinterp(:, i) = interp1(veczmds, Hmodes(:, i), MMP.z);
    
    eval(['umds' num2str(i-1) ' = Hmdsinterp(:, i) * u_mdsAmp(i, :);']);
    eval(['vmds' num2str(i-1) ' = Hmdsinterp(:, i) * v_mdsAmp(i, :);']);
    
    disp(['DONE WITH MODE ' num2str(i) ''])
end