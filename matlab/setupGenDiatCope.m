%
% Setup with generalists and a number of copepods
%
function p = setupGenDiatCope(n,nCopepods,mAdult, bParallel)

arguments
    n = 10;
    nCopepods = 10;
    mAdult (1,:) = [];
    bParallel = false;
end

loadNUMmodelLibrary();
calllib(loadNUMmodelLibrary(), 'f_setupgendiatcope',...
    int32(n), int32(nCopepods), length(mAdult), mAdult );
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupgendiatcope',....
            int32(n), int32(nCopepods),length(mAdult), mAdult );
    end
end

p.idxN = 1;
p.idxDOC = 2;
p.idxSi = 3;
p.idxB = 4; % We have three nutrient groups so biomass groups starts at index 4.

p.n = 4;
% Generalists
p = parametersAddgroup(5,p,n);
% Diatoms:
p = parametersAddgroup(3,p,n);

for i = 1:length(mAdult)
    p = parametersAddgroup(10,p,nCopepods, mAdult(i));
end


%p.n = p.n+p.idxB-1;
%p.nGroups = 1+length(mAdult);
%p.typeGroups = ones(1,1+p.nGroups);
%p.typeGroups(2:end) = 10;
%p.ixStart = p.idxB;
%p.ixEnd = p.n;

p = getMass(p);

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations)
p.u0(p.idxB:p.n) = 1;