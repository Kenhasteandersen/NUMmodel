%
% Sets up a spectrum of generalists simple.
%
% In:
%  n - number of size classes
%  bParallel - Whether to prepare parallel execution (for global runs)
%
function p = setupGeneralistsSimpleOnly(n, bParallel)

arguments
    n int32 {mustBeInteger, mustBePositive} = 10;
    bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

calllib(loadNUMmodelLibrary(), 'f_setupgeneralistssimpleonly', int32(n) );
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupgeneralistssimpleonly',int32(n));
    end
end

% Nutrients:
p = setupNutrients_N_DOC;

% Generalists:
p = parametersAddgroup(1,p,n);

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions
p.u0(p.idxB:p.n) = 1;