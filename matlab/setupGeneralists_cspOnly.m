function p = setupGeneralists_cspOnly(bParallel)

arguments
    %n int32 {mustBeInteger, mustBePositive} = 10;
    bParallel logical = false;
end

n = 10;

loadNUMmodelLibrary(bParallel);

calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly_csp' );
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly_csp');
    end
end

% Nutrients:
p = setupNutrients_N_DOC;

% Generalists:
p = parametersAddgroup(2,p,n);

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions (and deep layer concentrations)
p.u0(p.idxB:p.n) = 1;

% This model is usually run with quadratic HTL mortality:
parametersHTL(0.003, true)