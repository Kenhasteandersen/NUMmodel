%
% Setup with downregulating generalists and diatoms,  and a number of copepods
%
function p = setupGenDiatCope(mAdult,n,nCopepods,nPOM, bParallel)

arguments
    mAdult (1,:) = [];
    n = 10;
    nCopepods = 10;
    nPOM = 10;
    bParallel = false;
end

errortext ='                    ';
errorio=false;

[~,errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgendiatcope', ...
    int32(n), int32(nCopepods), int32(nPOM),length(mAdult), mAdult, errorio, errortext);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='                    ';
        [~,errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly', ...
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

% Diatoms :
p = parametersAddgroup(3,p,n);
% Generalists :
p = parametersAddgroup(5,p,n);


for i = 1:length(mAdult)
    p = parametersAddgroup(10,p,nCopepods, mAdult(i));
end

% POM:
p = parametersAddgroup(100, p, nPOM);

p = getMass(p);

p.u0(1:3) = [150, 0,200]; % Initial conditions (and deep layer concentrations)
% Initial condition at a Sheldon spectrum of "0.1":
% ix = 4:p.n;
% p.u0(ix) = 0.1*log( p.mUpper(ix)./p.mLower(ix)); 

p.u0(p.ixStart(1):p.ixEnd(end)) = 1;
p.u0( p.ixStart(end):p.ixEnd(end) ) = 0;% No POM in initial conditions

