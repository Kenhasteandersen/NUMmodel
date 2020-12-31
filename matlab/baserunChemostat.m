p = parameters(100);
%setParameters(p);

p = parametersChemostat(p);
p.tEnd = 365;

tic
sim = simulateChemostat(p, 100);
toc

plotSimulation(sim)