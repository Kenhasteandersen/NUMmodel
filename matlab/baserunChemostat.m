%
% Make a basic run of the chemostat model
%
p = parameters(0.1);

p = parametersChemostat(p);
p.tEnd = 365;

tic
p.bUseLibrary = false;
sim = simulateChemostat(p, 100);
toc

plotSimulation(sim)