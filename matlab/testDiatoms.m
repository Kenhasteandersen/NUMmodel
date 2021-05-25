function sim = testDiatoms()

%
% Set parameters:
%
p = setupDiatomsOnly;
p = parametersChemostat(p);
p.tEnd = 10;
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
sim = simulateChemostat(p, 60);
toc
%
% Plot
%
%plotChemostat(sim)