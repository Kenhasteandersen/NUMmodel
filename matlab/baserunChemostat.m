p = parameters(25);
setParameters(p);

p = parametersChemostat(p);
tic
sim = simulateChemostat(p, 10, 60);
toc
clf
loglog(p.m, sim.B(end,:))
ylim([1 1000])
