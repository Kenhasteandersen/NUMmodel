%
% Run a watercolumn with only generalists
%
function sim = baserunWatercolumn(lat,lon)

arguments
    lat double = 60;
    lon double = -10;
end

p = setupGeneralistsOnly(25);
p = parametersWatercolumn(p);

sim = simulateWatercolumn(p, lat,lon);
%%

day = 170;

figure(1)
plotGlobalWatercolumnTime(sim);

figure(2)
plotGlobalWatercolumn(sim,day);

figure(3)
plotGlobalSizespectrum(sim,day,1);

figure(4)
plotGlobalSizespectrumTime(sim,1);