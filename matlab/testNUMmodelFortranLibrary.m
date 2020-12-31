unloadlibrary(loadNUMmodelLibrary());
loadNUMmodelLibrary();

u(1) = 150.;
u(2) = 1.;
u(3:12) = 1:10;
dudt = 0*u;
L = 100;
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', 12, u, L, 0.1, dudt);


p = parameters(1e9);
rates = calcDerivatives(p, u, L);

format longe
rates.F./p.m

rates.dudt
dudt

%u = [p.N0 p.DOC0 ones(1,p.n)];
%dudt = 0*u;
%[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcrates',10., 60.,p.n+2,u, 1.0, 1.0, dudt);
%dudt


%p = parametersChemostat(p);
%tic
%sim = simulateChemostat(p, 10, 60);
%toc
%clf
%loglog(p.m, sim.B(end,:))
%ylim([1 1000])
