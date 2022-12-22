lat = 60;
lon = -10;
mAdult = logspace(log10(0.2), log10(1000), 7);
n = 10;
nCopepods = 10;
nPOM = 10;
p = setupGenDiatCope(mAdult, n,nCopepods,nPOM);
p = parametersWatercolumn(p);
p.tEnd = 365 * 10;
setHTL(0.15, 1, true, true);
sim = simulateWatercolumn(p, lat,lon); 