function p = setupGeneralistsDiatoms_simple(n, bParallel, options)

arguments
   n int32 {mustBeInteger, mustBePositive} = 10; % Number of grid points
   bParallel logical = false;
   options.bTest logical = false;
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

% Nutrients:
p = setupNutrients_N_DOC_Si;

% Generalists
p = parametersAddgroup(1,p,n);
% Diatoms simple:
p = parametersAddgroup(4,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 10]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart(1):p.ixEnd(end)) = 1;

if options.bTest
    u = p.u0;
    dudt = 0*u';
    [u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
        u, 60, 15, 0.1, dudt);
    if round(sum(dudt)*1000) == 1915
        p.bTestSuccess = true;
    else
        p.bTestSuccess = false;
    end
    
end