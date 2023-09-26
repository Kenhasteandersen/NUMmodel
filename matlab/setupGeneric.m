%
% Setup with generalists and a number of copepods
%
function p = setupGeneric(mAdult, bParallel)

arguments
    mAdult (1,:) = [];
    bParallel = false;
end

loadNUMmodelLibrary(bParallel);

errortext ='                    ';
errorio=false;

[~,errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult,errorio, errortext);
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='                    ';
        [~,errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult, errorio(i), this_errortext);
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

p = setupNutrients_N_DOC;

p.n = 2;
% Generalists:
p = parametersAddgroup(1,p,10);

for i = 1:length(mAdult)
    p = parametersAddgroup(10,p,10, mAdult(i));
end


%p.n = p.n+p.idxB-1;
%p.nGroups = 1+length(mAdult);
%p.typeGroups = ones(1,1+p.nGroups);
%p.typeGroups(2:end) = 10;
%p.ixStart = p.idxB;
%p.ixEnd = p.n;

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions
p.u0(p.idxB:p.n) = 1;
