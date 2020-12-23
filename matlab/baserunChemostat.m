p = parameters(10);
%setParameters(p);

p = parametersChemostat(p);
p.tEnd = 500;

tic
sim = simulateChemostat(p, 150);
toc

plotSimulation(sim)