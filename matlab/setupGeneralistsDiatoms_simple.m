function p = setupGeneralistsDiatoms_simple(n, bParallel)

arguments
   n int32 {mustBeInteger, mustBePositive} = 10; % Number of grid points
   bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsdiatoms_simple',int32(n));
    end
    p.bParallel = true;
else
    p.bParallel = false;
    calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsdiatoms_simple', int32(n) );
end

p.idxN = 1;
p.idxDOC = 2;
p.idxSi = 3;
p.idxB = 4; % We have three nutrient groups so biomass groups starts at index 4.

p.n = 3;
% Generalists
p = parametersAddgroup(1,p,n);
% Diatoms simple:
p = parametersAddgroup(4,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart(1):p.ixEnd(end)) = 1;