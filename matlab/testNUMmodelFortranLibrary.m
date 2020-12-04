p = parameters(10);
%%
unloadlibrary(loadNUMmodelLibrary());
setParameters(p);


u = [p.N0 p.DOC0 ones(1,p.n)];
dudt = 0*u;
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcrates',10., 60.,p.n+2,u, 1.0, 1.0, dudt);
dudt


%%
p = parametersChemostat(p);
tic
sim = simulateChemostat(p, 10, 60);
toc
clf
loglog(p.m, sim.B(end,:))
ylim([1 1000])
