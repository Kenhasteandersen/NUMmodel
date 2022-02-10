%
% Make a basic run of the chemostat model
%
function sim = baserunChemostatEuler(mAdult)
arguments
    mAdult double = []
end
%
% Set parameters:
%
p = setupGeneric(mAdult);
p = parametersChemostat(p);
p.tEnd = 365;
%
% Set to "normal" HTL mortality if there are no copepods:
%
if isempty(mAdult)
    setHTL(0.1, 1/500^1.5,false,false);
end
%
% Simulate
%
tic
sim = simulateChemostatEuler(p, 100);
toc
%
% Plot
%
plotSimulation(sim)