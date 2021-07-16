function sim = testGeneralistsDiatoms()

%%
% Set parameters:
%
p = setupGeneralistsDiatoms_simple;
p = parametersChemostat(p);

%
% Calc derivatives:
%
dudt = 0*ones(1,p.n);
u = p.u0;
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
    length(u), u, 60, 10., 0.0, dudt);
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
figure(1)
plotChemostat(sim)

%%
% Global test:
%
p = parametersGlobal(setupGeneralistsDiatoms_simple(10,true));
sim = simulateGlobal(p);

figure(2)
plotGlobal(sim)

figure(3)
plotWatercolumn(sim, 60, -10, 180)
