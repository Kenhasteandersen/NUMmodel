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
    lon double = -15;
end

close all force;
clc;
p = setupGeneralistsDiatoms();

p = parametersWatercolumn(p);
p.tEnd = 365;

sim = simulateWatercolumn(p, lat,lon);

exploreSimulations2(sim);

%plotSimulationTest(sim)

%checkConservation(sim);
