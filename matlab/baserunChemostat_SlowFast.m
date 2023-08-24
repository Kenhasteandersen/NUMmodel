% θ
% Make a basic run of the chemostat model
% In:
%  mAdult is the adult sizes of copepods (can be left empty to simulate only
%         unicellular organisms) (default = [], ie only generalists).
%
% Out:
%  sim: Structure holding the results of the simulation
%
function sim = baserunChemostat_SlowFast(mAdult)

arguments
    mAdult double = []
end

%Set parameters:

n = 10;
k = 2;
bParallel = false;

p = setupGeneralistsSimpleK(n,k, bParallel);

p = parametersChemostat(p);
p.tEnd = 200;
p.d = 0.1;

% Change colors for groups of Generalists:
j=linspace(0.3,0.9,k);
for i=1:k
    p.colGroup{i}=[0,0,j(i)];
end
%Change color of N so I can see generalists

p.colNutrients{1}=[0,1,1];

%p.umin=p.umin.*0.5;
%p.u0(p.idxB:p.n) = 0.5;
%
% Set to "normal" HTL mortality if there are no copepods:
%
if isempty(mAdult)
    setHTL(0.1, 1/500^1.5, false, false);
else 
    setHTL(0.1, 1, true, true);
end
%
% Simulate
%
tic
sim = simulateChemostat(p, 100);
toc
%
% Plot
%
plotSimulation_SlowFast(sim);
checkConservation(sim);
