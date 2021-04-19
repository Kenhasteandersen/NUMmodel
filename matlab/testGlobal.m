function bSuccess = testGlobal

p = parametersGlobal(setupGeneralistsOnly); % Use standard low-res model
p.tEnd = 5;
p.tSave = 5;

sim = simulateGlobal(p);

sumB = sum(sim.B(~isnan(sim.B)));
if ( sumB > 2.8e6 || sumB < 2.7e6 )
    bSuccess = false;
else
    bSuccess = true;
end

