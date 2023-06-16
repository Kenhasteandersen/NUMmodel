%
% Setup with generalists and a number of copepods
%
function p = setupNUMmodelSimple(mAdult, n,nCopepods,nPOM, options)

arguments
    mAdult (1,:) = [0.1 1 10 100 1000];
    n = 10;
    nCopepods = 10;
    nPOM = 10;
    options.bParallel = false;
end

loadNUMmodelLibrary(options.bParallel);

errortext ='';
errorio=false;

[~,errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupnummodelsimple', ...
    int32(n), int32(nCopepods), int32(nPOM),length(mAdult), mAdult,errorio, errortext);
if options.bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='';
        [~,errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupnummodelsimple', ...
            int32(n), int32(nCopepods), int32(nPOM),length(mAdult), mAdult, errorio(i), this_errortext);
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



p = setupNutrients_N_DOC_Si;
% Generalists simple:
p = parametersAddgroup(1,p,n);
% Diatoms simple:
p = parametersAddgroup(4,p,n);

for i = 1:length(mAdult)
    p = parametersAddgroup(11,p,nCopepods, mAdult(i)); % Active copepods
end

% POM:
p = parametersAddgroup(100, p, nPOM);

p = getMass(p);

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations)
% Initial condition at a Sheldon spectrum of "0.1":
ix = p.idxB:p.n;
p.u0(ix) = 0.1*log( p.mUpper(ix)./p.mLower(ix)); 

p.u0( p.ixStart(end):p.ixEnd(end) ) = 0; % No POM in initial conditions