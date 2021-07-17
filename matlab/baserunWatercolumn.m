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

p = setupGeneralistsOnly(25);
p = parametersWatercolumn(p);

sim = simulateWatercolumn(p, lat,lon);
%%

day = 170;

figure(1)
plotWatercolumnTime(sim,'depthMax',200);

figure(2)
plotWatercolumn(sim,day,'depthMax',200);

figure(3)
plotGlobalSizespectrum(sim,day,1);

figure(4)
plotSizespectrumTime(sim,1);