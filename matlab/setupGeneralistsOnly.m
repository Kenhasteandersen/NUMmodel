%
% Sets up a spectrum of generalists.
%
% In:
%  n - number of size classes
%  bParallel - Whether to prepare parallel execution (for global runs)
%
function p = setupGeneralistsOnly(n, bParallel)

arguments
    n int32 {mustBeInteger, mustBePositive} = 10;
    bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

errortext ='                    ';
errorio=false;
[errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly', int32(n), errorio, errortext );

if errorio
    disp(['Error loading ',errortext,'. Execution terminated'])
    return
else
    disp('done loading input parameters')
end


if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;%f_setupgeneralistssonly

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);
    parfor i=1:poolsize
        this_errortext ='                    ';
        [errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly', int32(n), errorio(i), this_errortext);
        errortext(i)={this_errortext}
    end
    if any(errorio)
        i=find(errorio==true,1);
        disp(['Error loading ',errortext{i},'. Execution terminated'])
        return
    else
        disp('done loading input parameters')
    end
end


% Nutrients:
p = setupNutrients_N_DOC;

% Generalists:
p = parametersAddgroup(5,p,n);

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions
p.u0(p.idxB:p.n) = 1;
