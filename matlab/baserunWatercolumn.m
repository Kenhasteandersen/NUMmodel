%
% Run a watercolumn with only generalists
%
% In:
%  lat, lon - latitude and longitude
%
% Out:
%  As simulation structure
%
function sim = baserunWatercolumn(lat, lon)

arguments
    lat double = 0;
    lon double = -172;
end

p = setupNUMmodel();

p = parametersWatercolumn(p);
p.tEnd = 5*365;

setHTL(0.1, 1, true, true)
setSinkingPOM(p, 0.78622);
sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);
