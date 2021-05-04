function p = setupGeneralistsOnly(n, bParallel)

arguments
    n int32 {mustBeInteger, mustBePositive} = 10;
    bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        %loadNUMmodelLibrary();
        calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly',int32(10));
    end
else
    calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly', int32(n) );
end

p.idxN = 1;
p.idxDOC = 2;
p.idxB = 3; % We have two nutrient groups so biomass groups starts at index 3.

p.n = 2;
% Generalists:
p = parametersAddgroup(1,p,n);

p = getMass(p);
%[p.m(p.idxB:p.n), p.mLower(p.idxB:p.n), p.mDelta(p.idxB:p.n)] = parametersCalcGrid(10^-8.5, 0.1, n);

p.u0(1:2) = [150, 0]; % Initial conditions (and deep layer concentrations)
p.u0(p.idxB:p.n) = 1;