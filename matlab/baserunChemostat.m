%
% Make a basic run of the chemostat model
% In:
%  mAdult is the adult sizes of copepods (can be left empty to simulate only
%         unicellular organisms) (default = [], ie only generalists).
%  bUseFortran: Flag indicating whether the fortran library is used
%  (default to false)
%
% Out:
%  sim: Structure holding the results of the simulation
%
function sim = baserunChemostat(mAdult)

arguments
    mAdult double = []
end
    
%
% Set parameters:
%
p = setupGeneralistsOnly();
p = parametersChemostat(p);
p.tEnd = 365;
%
% Simulate
%
tic
sim = simulateChemostat(p, 100);
toc
%
% Plot
%
plotChemostat(sim)