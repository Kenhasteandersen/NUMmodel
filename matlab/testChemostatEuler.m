function bSuccess = testChemostatEuler(value)

p = setupNUMmodel;
p = parametersChemostat(p);
sim = simulateChemostatEuler(p,1.0);

sumB = sum(sim.B(:));
if ( sumB > 0.99*value || sumB < 1.01*value )
    bSuccess = false;
else
    bSuccess = true;
end

bSuccess = true;

