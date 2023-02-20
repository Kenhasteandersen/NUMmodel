%
% Setup with generalists and a number of copepods
%
function p = setupNUMmodel(mAdultPassive, mAdultActive, n,nCopepods,nPOM, bParallel)

arguments
    mAdultPassive (1,:) = [0.2 5];
    mAdultActive (1,:) = [1 10 100 1000];
    n = 10;
    nCopepods = 10;
    nPOM = 10;
    bParallel = false;
end

loadNUMmodelLibrary(bParallel);
calllib(loadNUMmodelLibrary(), 'f_setupnummodel', ...
    int32(n), int32(nCopepods), int32(nPOM), ...
    length(mAdultPassive), mAdultPassive, length(mAdultActive), mAdultActive );
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupnummodel', ...
            int32(n), int32(nCopepods), int32(nPOM),...
            length(mAdultPassive), mAdultPassive, length(mAdultActive), mAdultActive );
    end
end

p = setupNutrients_N_DOC_Si;

% Generalists:
p = parametersAddgroup(5,p,n);

% Diatoms
p = parametersAddgroup(3,p,n);

% Copepods:
for i = 1:length(mAdultPassive)
    p = parametersAddgroup(10,p, nCopepods, mAdultPassive(i));
end
for i = 1:length(mAdultActive)
    p = parametersAddgroup(11,p, nCopepods, mAdultActive(i));
end

% POM:
p = parametersAddgroup(100, p, nPOM);

p = getMass(p);

p.u0(1:3) = [150, 0, 10]; % Initial conditions (and deep layer concentrations)
% Initial condition at a Sheldon spectrum of "0.1":
ix = p.idxB:p.n;
p.u0(ix) = 0.1*log( p.mUpper(ix)./p.mLower(ix)); 

p.u0( p.ixStart(end):p.ixEnd(end) ) = 0; % No POM in initial conditions