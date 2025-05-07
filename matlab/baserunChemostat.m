% 
% Make a basic run of the chemostat model using the NUMmodel setup
%
% Out:
%  sim: Structure holding the results of the simulation
%


p = setupGeneralistsOnly();       % Sets up the model
p = parametersChemostat(p);% Sets up the chemostat environment

p.tEnd = 500; % Run 500 days
p.d = 0.05;   % Set the mixing rate to 0.05/day
%
% Simulate
%
tic
sim = simulateChemostat(p, 50);
toc
%
% Plot
%
plotSimulation(sim);
checkConservation(sim);
