%
% Setup with generalists, diatoms, passive a number of specific copepod species, and POM
%
% In:
%  - idSpecies: an array of the species numbers. These numbers should
%               correspond to species numbers in input.h
%
function p = setupNUMmodelSpecies(idSpecies, n,nCopepods,nPOM, options)

arguments
    idSpecies = 1;
    n = 10;  % Number of size groups in generalist and diatom spectra
    nCopepods = 10;  % Number of stages in copepod groups
    nPOM = 1;  % Number of POM size groups
    options.bParallel = false;  % Whether to prepare for parallel runs (for global runs)
end
%
% Load the NUM library and call the setupNUMmodel:
% 
loadNUMmodelLibrary(options.bParallel);

errortext ='                    ';
errorio=false;


%calllib(loadNUMmodelLibrary(), 'f_testsetup', 11.4)

[~,errorio,errortext]= ...
calllib(loadNUMmodelLibrary(), 'f_setupnummodelspecies',...
    int32(n), int32(nCopepods), int32(nPOM),...
    int32(length(idSpecies)), int32(idSpecies), ...
    errorio, errortext);
if options.bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='                    ';
        [errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupnummodelspecies', ...
            int32(n), int32(nCopepods), int32(nPOM),...
            length(idSpecies), int32(idSpecies), ...
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
for i = 1:length(idSpecies)
    p = parametersAddgroup(49+idSpecies(i), p, nCopepods, 0);
end

% POM:
p = parametersAddgroup(100, p, nPOM);
setSinkingPOM(p, 20); 

% Initial conditions:
p = getMass(p);

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations) of nutrients
% Initial condition at a Sheldon spectrum of "0.1":
ix = p.idxB:p.n;
p.u0(ix) = 0.1*log( p.mUpper(ix)./p.mLower(ix)); 

p.u0( p.ixStart(end):p.ixEnd(end) ) = 0; % No POM in initial conditions

setHTL(0.015,0.1,true,false);
