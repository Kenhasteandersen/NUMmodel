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
    lat double = 43;%60 %43;49 ; 31
    lon double = 8;%-40 %8;-4 ; -64
end

p = setupNUMmodelzoo();

p = parametersWatercolumn(p);
p.tEnd = 5*365;

sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);
end
