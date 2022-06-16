% 
% Make a basic run of the new implementation of POM matters - Author:
% Cécile Decker
% 
% In:
%
% Out:
%  sim: Structure holding the results of the simulation
%
function sim = baserunPOM()

mAdult = logspace(log10(0.2), log10(1000), 5);
n = 10;
nCopepods = 10;
nPOM = 10;
p = setupNUMmodel(mAdult, n,nCopepods,nPOM);

p = parametersWatercolumn(p);
p.tEnd = 365*2;

lat = 60;
lon = -10;
sim = simulateWatercolumn(p, lat,lon);

if strcmp(sim.p.nameModel,'watercolumn') || strcmp(sim.p.nameModel,'global')
    sim = calcFunctions(sim);
end

sim = calcPOM(sim);

plotPOM(sim);