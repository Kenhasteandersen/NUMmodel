function p = setupDiatoms_simpleOnly(n, bParallel)

arguments
   n int32 {mustBeInteger, mustBePositive} = 10; % Number of grid points
   bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);
    parfor i=1:poolsize
        this_errortext ='                    ';
        [errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupdiatoms_simpleonly',int32(n),errorio(i), this_errortext);
        errortext(i)={this_errortext}
    end
    p.bParallel = true;
    if any(errorio)
        i=find(errorio==true,1);
        disp(['Error loading ',errortext{i},'. Execution terminated'])
        return
    else
        disp('done loading input parameters')
    end
else
    errortext ='                    ';
    errorio=false;
    [errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupdiatoms_simpleonly', int32(n),errorio, errortext);
    p.bParallel = false;
    if errorio
        disp(['Error loading ',errortext,'. Execution terminated'])
        return
    else
        disp('done loading input parameters')
    end
end

% Nutrients:
p = setupNutrients_N_DOC_Si;

% Diatoms:
p = parametersAddgroup(4,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 10]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart:p.ixEnd) = 1;
