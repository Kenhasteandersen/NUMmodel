function p = setupDiatomsOnly(n, bParallel)

arguments
   n int32 {mustBeInteger, mustBePositive} = 10; % Number of grid points
   bParallel logical = false;
end

loadNUMmodelLibrary(bParallel);

errortext ='                    ';
errorio=false;

[errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly', int32(n),errorio, errortext);

if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;

    errorio=false(1,poolsize);
    errortext = repmat({''}, [1 poolsize]);

    parfor i=1:poolsize
        this_errortext ='                    ';
        [errorio(i),this_errortext]=calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly',int32(n),errorio(i), this_errortext);
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


%**************** previous version ********************
%
% if bParallel
%     h = gcp('nocreate');
%     poolsize = h.NumWorkers;
%     parfor i=1:poolsize
%         calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly',int32(n));
%     end
%     p.bParallel = true;
% else
%     calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly', int32(n) );
%     p.bParallel = false;
% end
%
%*************** end of prev version **************

% Nutrients:
p = setupNutrients_N_DOC_Si;

% Diatoms:
p = parametersAddgroup(3,p,n);

p = getMass(p); % Get masses

p.u0(1:3) = [150, 0, 200]; % Initial conditions (and deep layer concentrations)
p.u0(p.ixStart:p.ixEnd) = 1;
