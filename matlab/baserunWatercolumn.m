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

%     lat double = 0;
%     lon double = -23;

end

mAdults = logspace(log10(0.2), log10(1000), 10);


% p = setupGeneralistsOnly(25);
p = setupGeneric(mAdults);
p = parametersWatercolumn(p);
% p.tEnd = 2*365;
p.tEnd = 1000;

setHTL(0.1, 10, true, true);

sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);