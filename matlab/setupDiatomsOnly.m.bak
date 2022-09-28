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

% Nutrients:
p = setupNutrients_N_DOC_Si;

% Diatoms:
p = parametersAddgroup(3,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 10]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart:p.ixEnd) = 1;