%
% Run a watercolumn with only generalists
%
% In:
%  lat, lon - latitude and longitude
%
% Out:
%  As simulation structure
%
function sim = baserunWatercolumn(mAdult, lat, lon)

arguments
    mAdult double = []
    lat double = 60;
    lon double = -10;
end

% p = setupGeneralistsOnly(25);
p = setupGeneric(mAdult);
p = parametersWatercolumn(p);
p.tEnd = 1095;

%
% Set to "normal" HTL mortality if there are no copepods:
%
if isempty(mAdult)
    setHTL(0.1, 1/500^1.5, false, false);
else 
    setHTL(0.1, 1, true, true);
end

sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);