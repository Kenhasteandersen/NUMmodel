%
% Setup with generalists, diatoms, passive and active copepods, and POM.
%
% In:
%   mAdultPassive - the adult masses of passive copepods (default [0.2 5])
%   mAdultActive  - the adult masses of active copepods (logspace(0,3,3) )
%   n             - No. of unicellular size groups (10)
%   nCopepods     - No. of stages in the copepod groups (6)
%   nPOM          - No. of POM size groups (1)
%   bParallel     - whether to prepare a parallel run (FALSE)
%
% Out:
%   p             - The parameter structure
%
function p = setupNUMmodelzoo(mAdultGelatinous, mAdultNongelatinous, n, nZooplanktons, nPOM, options)

arguments
    mAdultGelatinous (1,:) = logspace(0,6,6) %  [0.2 5];  % Adult masses of passive copepods
    mAdultNongelatinous (1,:) = logspace(0,6,6);  % 3 log-spaced adult masses of active copepods
    n = 10;  % Number of size groups in generalist and diatom spectra
    nZooplanktons = 6;  % Number of stages in copepod groups
    nPOM = 1;  % Number of POM size groups
    options.bParallel = false;  % Whether to prepare for parallel runs (for global runs)
end
%
% Load the NUM library and call the setupNUMmodel:
% 
loadNUMmodelLibrary(options.bParallel);

errortext ='                    ';
errorio=false;

[~,~,errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupnummodelzoo', ...
    int32(n), int32(nZooplanktons), int32(nPOM), ...
    length(mAdultGelatinous), mAdultGelatinous, length(mAdultNongelatinous), mAdultNongelatinous, ...
    errorio, errortext);
if options.bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='                    ';
        [~,~,errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupnummodelzoo', ...
            int32(n), int32(nZooplanktons), int32(nPOM),...
            length(mAdultGelatinous), mAdultGelatinous, length(mAdultNongelatinous), mAdultNongelatinous ,...
            errorio(i), this_errortext);
        errortext(i)={this_errortext}
    end
    if any(errorio)
        i=find(errorio==true,1);
        disp(['Error loading ',errortext{i},'. Execution terminated'])
        return
    %else
    %    disp('done loading input parameters')
    end
else
    if errorio
        disp(['Error loading ',errortext,'. Execution terminated'])
        return
    %else
    %    disp('done loading input parameters')
    end
end
%
% Setup in matlab:
%
p = setupNutrients_N_DOC_Si;

% Generalists:
p = parametersAddgroup(5,p,n);

% Diatoms
p = parametersAddgroup(3,p,n);

% Copepods:
for i = 1:length(mAdultGelatinous)
    p = parametersAddgroup(12,p, nZooplanktons, mAdultGelatinous(i));
end
for i = 1:length(mAdultNongelatinous)
    p = parametersAddgroup(13,p, nZooplanktons, mAdultNongelatinous(i));
end

% POM:
p = parametersAddgroup(100, p, nPOM);
setSinkingPOM(p, 19); 

% Initial conditions:
p = getMass(p);

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations) of nutrients
% Initial condition at a Sheldon spectrum of "0.1":
ix = p.idxB:p.n;
p.u0(ix) = 0.1*log( p.mUpper(ix)./p.mLower(ix) ); 

p.u0( p.ixStart(end):p.ixEnd(end) ) = 0; % No POM in initial conditions

setHTL(0.017, 1 ,true, false, false,true); % "Quadratic" mortality; not declining; only affecting copepods
