function [MdsProf, MdsAmp, VarMds] = modalDecomp(nmds, zN2, N2, D, varStruct, N2tandwnd)
% [MdsProf, MdsAmp, VarMds] = MODALDECOMPOSITION(zN2, N2, D, zH, varH, zV, var, N2window)
%
%   inputs:
%       - nmds: number of BAROCLINIC(???) modes to fit (FIX THIS PURPOSE!!!)
%       - zN2: as required by inertGravVmodes function.
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
% If N2 is a matrix and no window is specified, then normal modes
% are computed for every column of N2. This greatly increases the
% computational time.
%
% Olavo Badaro Marques, 18/Nov/2016.

tic

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

% Cell array with the name of all variables to decompose
% into modes and how many of them there are:
namesallvar = [namesVarH; namesVarV];
nvars = length(namesallvar);

lforHvar = [true(length(namesVarH), 1); false(length(namesVarV), 1)];
typesallvar = cell(1, nvars);

typesallvar(lforHvar) = {'varH'};   % string to indicate appropriate
typesallvar(~lforHvar) = {'varV'};  % structure field in a dynamical "way"
                                    % (get correct terminology, ???)

% Number of profiles we have to fit:
nprofs = size(varStruct.(typesallvar{1}).(namesallvar{1}), 2);
% FOR NOW (AND PROBABLY FOREVER), LET'S TAKE ALL THE DATA GRIDDED IN TIME.


%% Create/replicate z vector cell array for each variable:

% -------------------------------------------------------
% check that varStruct.z is not specified along with
% varStruct.H/varStruct.V, throw error if it is the case
%
% test to see if it works......

if isfield(varStruct, 'z') && ...
       (isfield(varStruct, 'zH') || isfield(varStruct, 'zV')) % parenthesis
                                                             % not required 
    error('NOT ALLOWED!!!')
    
else
    
    if ~isfield(varStruct, 'z')
        
        if ~isfield(varStruct, 'zH')
            varStruct.zH = {};
        end
        
        if ~isfield(varStruct, 'zV')
            varStruct.zV = {};
        end
        
    end
    
end


% Now create the zvecs cell-array, with the depth vectors associated
% with each of the variables to do the fit:

% Pre-allocate:
zvecs = cell(1, nvars);


if isfield(varStruct, 'z')    
    % Replicate varStruct.z to all fields of zvecs
    zvecs(:) = {varStruct.z};    
    
else
    
    nvarsH = length(namesVarH);
    nvarsV = length(namesVarV);
    
    if length(varStruct.zH)==1
        zvecs(1:nvarsH) = {varStruct.zH};
    else
        zvecs(1:nvarsH) = varStruct.zH;
    end
      
    
    if length(varStruct.zV)==1        
        zvecs(nvarsH+1:nvarsV) = {varStruct.zV};
    else
        zvecs(nvarsH+1:nvarsV) = varStruct.zH;
    end
    
end

% -------------------------------------------------------



%% Check whether N2tandwnd was specified and determine
% how many chuncks of data will be used:

% % % N2time = N2tandwnd{1};
% % % N2wndw = N2tandwnd{2};
% % % 
% % % N2halfwndw = round(N2wndw/2);

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
% THE WINDOWED CASE IS MISSING FROM THE IF ABOVE!!! WHICH SHOULD PROBABLY
% CREATE A DIFFERENT N2 MATRIX, IF I'M TO STICK WITH SUBSETTING WITH THE
% indN2col VARIABLE....


%% Pre-allocate space for the output variables:

% Loop through variable names:
for i1 = 1:length(namesallvar)
    
    % Check whether to include the "barotropic mode":
    if lforHvar(i1)
        vecmds = 0:nmds;
    else
        vecmds = 1:nmds;
    end
    
    % Pre-allocate space for modal structures output:
%     MdsProf. = NaN(length)
% I would love to compute the modal variables by multiplying structure and
% amplitude....except that different variables may have different z grid...
% Make such that the output is the zN2aug grid. However, it may even be
% faster to interpolate the modal structure at the end.
    
    % Pre-allocate space for modal amplitude output:
    MdsAmp.(namesallvar{i1}) = NaN(length(vecmds), nprofs);
    
    % Loop through modes and pre-allocate space for modal variables:
    for i2 = 1:length(vecmds)
        % NaN matrix with the same size as correspondent input:
        VarMds.(namesallvar{i1}).(['md' num2str(vecmds(i2)) '']) = ...
            NaN(size(varStruct.(typesallvar{i1}).(namesallvar{i1})));
    end
end


%% Computing modal amplitudes -- this is the hard part of the code:

zN2aug = [0; zN2; D];

% -------------------------------------------------------------------------
% If I'm going to loop through every profile anyway, I should definitely
% rearrange this loop. Well....if there are no gaps, I definitely don't
% want to loop through profiles:

if ndiffMds==nprofs
    % BRUTE FORCE - FIT FOR EVERY PROFILE
    
    for i1 = 1:nprofs
    
        [Hmodes, Vmodes, ~] = inertGravVmodes(zN2, D, ...
                                              N2(:, indN2col(i1)), nmds);

        for i2 = 1:nvars
            
            % Check which mode do we need to fit and use "Hmodes" or
            % "Vmodes(:, 2:end)". I can probably do something better than
            % this if....
            if lforHvar(i2)
                
                MdsAmp.(namesallvar{i2})(:, i1) = fitVmodes(zvecs{i2},   ...
                  varStruct.(typesallvar{i2}).(namesallvar{i2})(:, i1), ...
                  Hmodes, zN2aug);
              
            else
               
                MdsAmp.(namesallvar{i2})(:, i1) = fitVmodes(zvecs{i2},   ...
                  varStruct.(typesallvar{i2}).(namesallvar{i2})(:, i1), ...
                  Vmodes(:, 2:end), zN2aug);
                
            end
            
            
        end

        % Print progress to the screen:
        if mod(i1, round(0.1*nprofs))==0
            disp(['Fitting modes to column ' num2str(i1) ' of ' ...
                  '' num2str(nprofs) ''])
        end
    
      
    end
  
    
%% 
else


    for i1 = 1:ndiffMds
            
        % I should first subset the appropriate window (window case
        % is not implemented yet) FOR THE DATA!!!:
        
        % Find modes:
        [Hmodes, Vmodes, ~] = inertGravVmodes(zN2, D, ...
                                              N2(:, indN2col(i1)), nmds);
        zmds = [0; zN2; D];
                                          
        % Now I find which columns have the same NaNs:
        
        for i2 = 1:nvars
            
            if lforHvar(i2)
                vecmds = 0:nmds;
            else
                vecmds = 1:nmds;
            end
            
            % Must implement the efficient version later... or I could
            % simply loop through all columns and use fitVmodes
            [colsaux, rnanaux] = sepColsNanPos(varStruct.(typesallvar{i2}).(namesallvar{i2}));
            
            % SHOULD INCLUDE AN OPTION TO AVOID DOING THE ABOVE IN THE CASE
            % OF SAME GRID (THE ABOVE IS REDUNDANT FOR U AND V)
            
            % I have to interpolate modes for the appropriate zvecs{i2}
            vmodesinterp = NaN(length(zvecs{i2}), nmds);

            for imds = 1:(nmds+1)
                vmodesinterp(:, imds) = interp1(zmds, Hmodes(:, imds), ...
                                                    zvecs{i2});
            end
                
            for i3 = 1:length(colsaux)
                
                %
                rokgood = setdiff(1:length(zvecs{i2}), rnanaux{i3});
                % I could also ouput rokgood when outputting rnanaux
                
                % Subset modes for good data:
                vmds_aux = vmodesinterp(rokgood, :);
                
                m = ( vmds_aux' * vmds_aux) \ ( vmds_aux' * ...
                        varStruct.(typesallvar{i2}).(namesallvar{i2})(rokgood, colsaux{i3}) );
                
                MdsAmp.(namesallvar{i2})(:, colsaux{i3}) = m;
                
                for i4 = 1:length(vecmds)
                
                    VarMds.(namesallvar{i1}).(['md' num2str(vecmds(i4)) ''])(:, colsaux{i3}) = ...
                        vmodesinterp(:, i4) * m(i4, :);
                end
                
            end
                        
        end
                
    end
end
toc
 
% WHY IS THERE A COLUMN OF NAN for u somewhere...?????? (it is about at the
% maximum of mode-1, somewhere,  around yday 38)

% -------------------------------------------------------------------------

keyboard
