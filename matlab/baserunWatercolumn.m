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
    lon double = -40;
end

% p = setupDiatomsOnly(10);
mAdult = logspace(log10(0.2), log10(1000), 3);

n = 5;
nCopepods = 10;
nPOM = 10;
p = setupGenDiatCope(mAdult, n,nCopepods,nPOM);
% p = setupNUMmodel(mAdult, n,nCopepods,nPOM);

% p = setupGeneralistsDiatoms(10);
% p = setupGeneric(mAdult);
p = parametersWatercolumn(p);
p.tEnd = 5*365;

%
% Set to "normal" HTL mortality if there are no copepods:
%
% if isempty(mAdult)
%     setHTL(0.1, 1/500^1.5, false, false);
% else
%     setHTL(0.1, 1, true, true);
% end
setHTL(0.15, 1, true, true);

sim = simulateWatercolumn(p, lat,lon);

plotSimulation(sim)

checkConservation(sim);