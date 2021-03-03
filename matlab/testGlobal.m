function bSuccess = testGlobal

p = parameters([]);
p = parametersGlobal(p); % Use standard low-res model
p.tEnd = 5;
p.tSave = 5;

sim = simulateGlobal(p);

sumB = sum(sim.B(~isnan(sim.B)));
if ( sumB > 2830000 | sumB < 2820000 )
    bSuccess = false;
else
    bSuccess = true;
end

