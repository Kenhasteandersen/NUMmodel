function bSuccess = testGlobal(value)

p = setupGeneralistsPOM(5,1, true); % Fast setup with POM
p = parametersGlobal(p); % Use standard low-res model
p.tEnd = 30;
p.tSave = 10;

sim = simulateGlobal(p);
checkConservation(sim);
%fprintf("N conservation: %f per year\n", (sim.Ntot(end)/sim.Ntot(1)-1)/(p.tEnd/365) )

sumB = sum(sim.B(~isnan(sim.B)));
if ( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end


