function bSuccess = testGlobal(value)

p = setupGeneralistsOnly(10, true);
p = parametersGlobal(p); % Use standard low-res model
p.tEnd = 5;
p.tSave = 5;

sim = simulateGlobal(p);

sumB = sum(sim.B(~isnan(sim.B)));
if ( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

