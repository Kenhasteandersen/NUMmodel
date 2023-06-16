options.bParallel = false;
n = 10;
%Clost=0;
errorio=false;

loadNUMmodelLibrary(options.bParallel);
errortext ='init';

[errorio,errortext]=calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly', int32(n), errorio, errortext );

%% generalists simple
[errorio,text]=calllib(loadNUMmodelLibrary(), 'f_setupnummodel2',int32(n), errorio, text );
if errorio
    disp(['Error loading ',text,'. Execution terminated'])
    return
end
% [errorstr]=calllib(loadNUMmodelLibrary(), 'f_setupnummodel2',int32(n),  errorstr );
% calllib(loadNUMmodelLibrary(), 'f_setupgeneralistssimpleonly', int32(n));
%% generalists
% calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly', int32(n));





%% setup and run

p = setupNutrients_N_DOC;

% Generalists:
p = parametersAddgroup(1,p,n);

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions
p.u0(p.idxB:p.n) = 1;


p = parametersChemostat(p);
p.tEnd = 200;
p.d = 0.1;
tic
sim = simulateChemostat(p, 100);
toc

plotSimulation(sim);