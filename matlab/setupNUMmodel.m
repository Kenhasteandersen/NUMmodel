%
% Setup with generalists and a number of copepods
%
function p = setupNUMmodel(mAdult, n,nCopepods,nPOM, bParallel)

arguments
    mAdult (1,:) = [];
    n = 10;
    nCopepods = 10;
    nPOM = 10;
    bParallel = false;
end

loadNUMmodelLibrary();
calllib(loadNUMmodelLibrary(), 'f_setupnummodel', ...
    int32(n), int32(nCopepods), int32(nPOM),length(mAdult), mAdult );
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupnummodel', ...
            int32(n), int32(nCopepods), int32(nPOM),length(mAdult), mAdult );
    end
end

p.idxN = 1;
p.idxDOC = 2;
p.idxB = 3; % We have two nutrient groups so biomass groups starts at index 3.

p.n = 2;
% Generalists:
p = parametersAddgroup(1,p,n);

for i = 1:length(mAdult)
    p = parametersAddgroup(10,p,nCopepods, mAdult(i));
end

% POM:
p = parametersAddgroup(100, p, nPOM);

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions (and deep layer concentrations)
% Initial condition at a Sheldon spectrum of "0.1":
ix = 3:p.n;
p.u0(ix) = 0.1*log( p.mUpper(ix)./p.mLower(ix)); 

p.u0( p.ixStart(end):p.ixEnd(end) ) = 0; % No POM in initial conditions
