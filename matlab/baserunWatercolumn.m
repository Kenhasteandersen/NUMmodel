%
% Run a watercolumn with only generalists
%
% In:
%  lat, lon - latitude and longitude
%
% Out:
%  As simulation structure
%
function sim = baserunWatercolumn(lat,lon)

arguments
    lat double = 60;
    lon double = -10;
end

% p = setupGeneralistsOnly(25);
p = setupDiatoms_simpleOnly(10);
% p = setupDiatomsOnly;
p = parametersWatercolumn(p);
p.tEnd = 2*365;

sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);