%
% Run a watercolumn with the NUMmodel setup
%
% In:
%  lat, lon - latitude and longitude
%
% Out:
%  A simulation structure
%
function sim = baserunWatercolumn(lat, lon)

arguments
    lat double = 60;
    lon double = -40;
end

p = setupNUMmodel();

p = parametersWatercolumn(p);
p.tEnd = 5*365;

sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);
