%
% Sets up a spectrum of generalists with POM..
%
% In:
%  n - number of size classes for generalists
%  nPOM number of size classes for POM
%  bParallel - Whether to prepare parallel execution (for global runs)
%
function p = setupGeneralistsPOM(n, nPOM, bParallel)

arguments
    n int32 {mustBeInteger, mustBePositive} = 10;
    nPOM int32 {mustBeInteger, mustBePositive} = 10;
    bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

errortext ='                    ';
errorio=false;

[errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistspom', int32(n), int32(nPOM),errorio, errortext);


if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);
    parfor i=1:poolsize
        this_errortext ='                    ';
        [errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistspom', int32(n), int32(nPOM),errorio(i), this_errortext);
      errortext(i)={this_errortext}  
    end
    if any(errorio)
        i=find(errorio==true,1);
        disp(['Error loading ',errortext{i},'. Execution terminated'])
        return
    else
        disp('done loading input parameters')
    end
else
    if errorio
        disp(['Error loading ',errortext,'. Execution terminated'])
        return
    else
        disp('done loading input parameters')
    end
end


% Nutrients:
p = setupNutrients_N_DOC;

% Generalists:
p = parametersAddgroup(5,p,n);
p.u0(p.ixStart(1):p.ixEnd(1)) = 1; % Initial conditions

% POM:
p = parametersAddgroup(100, p, nPOM);
p.u0(p.ixStart(2):p.ixEnd(2)) = 0; % Initial conditions

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions (and deep layer concentrations)
