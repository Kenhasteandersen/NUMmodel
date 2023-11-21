%
% Sets up a spectrum of generalists simple.
%
% In:
%  n - number of size classes
%  bParallel - Whether to prepare parallel execution (for global runs)
%
function p = setupGeneralistsSimpleK(n,k, bParallel)

arguments
    n int32 {mustBeInteger, mustBePositive} = 10;
    k int32 {mustBeInteger, mustBePositive} = 2;
    bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);
errortext ='';
errorio=false;


[errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistssimple_two', int32(n), int32(k),errorio, errortext);
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='';
        [errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistssimple_two', int32(n), int32(k),errorio(i), this_errortext);
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
% Prokaryote
p = parametersAddgroup(6,p,n);

% Generalists:
for i=1:k
p = parametersAddgroup(1,p,n);
end

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions
p.u0(p.idxB:p.n) = 1;