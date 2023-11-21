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
k = 3;
bParallel = false;

p = setupGeneralistsSimpleK(n,k, bParallel);

p = parametersChemostat(p);
p.tEnd = 720;
p.d = 0.1;

% Change colors for groups of Generalists:
r=[13 15 56 119 189];
g=[139 178 192 207 227];
b=[217 242 242 242 242];
j=linspace(0.3,0.9,k);

blue(1,:)=[58/255,67/255,186/255];
blue(2,:)=[4/255,146/255,194/255];
blue(3,:)=[130/255,237/255,253/255];

for i=2:k+1
    %p.colGroup{i}=[0,0,j(i-1)];
    p.colGroup{i}=blue(i-1,:);
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
% plotSimulation_SlowFast(sim);
%plotSimulation(sim);
plotSimulation_trine(sim);
checkConservation(sim);
