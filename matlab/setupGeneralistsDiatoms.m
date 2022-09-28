function p = setupGeneralistsDiatoms(n, bParallel)

arguments
   n int32 {mustBeInteger, mustBePositive} = 10; % Number of grid points
   bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsdiatoms',int32(n));
    end
else
    calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsdiatoms', int32(n) );
end

% Nutrients:
p = setupNutrients_N_DOC_Si;

% Generalists
p = parametersAddgroup(5,p,n);
% Diatoms:
p = parametersAddgroup(3,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart(1):p.ixEnd(end)) = 1;