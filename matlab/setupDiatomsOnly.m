function p = setupDiatomsOnly(n, bParallel)

arguments
   n int32 {mustBeInteger, mustBePositive} = 10; % Number of grid points
   bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly',int32(n));
    end
    p.bParallel = true;
else
    calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly', int32(n) );
    p.bParallel = false;
end

p.idxN = 1;
p.idxDOC = 2;
p.idxSi = 3;
p.idxB = 4; % We have three nutrient groups so biomass groups starts at index 4.

p.n = 3;
% Diatoms:
p = parametersAddgroup(3,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart:p.ixEnd) = 1;