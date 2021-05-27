function sim = testGeneralistsDiatoms()

%%
% Set parameters:
%
p = setupGeneralistsDiatoms;
p = parametersChemostat(p);

%
% Calc derivatives:
%
dudt = 0*ones(1,p.n);
u = p.u0;
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
    length(u), u, 60, 0.0, dudt);
dudt
%
% Simulate
%
tic
p.tEnd = 200;
sim = simulateChemostat(p, 60);
toc
%
% Plot
%
plotChemostat(sim)

%%
% Global test:
%
%p = parametersGlobal(setupDiatomsOnly(10,false));
%p.tEnd = 10;
%sim = simulateGlobal(p);