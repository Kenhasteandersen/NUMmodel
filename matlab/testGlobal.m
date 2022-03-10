function bSuccess = testGlobal

p = setupGeneralistsOnly(10, true);
p = parametersGlobal(p); % Use standard low-res model
p.tEnd = 5;
p.tSave = 5;

sim = simulateGlobal(p);

sumB = sum(sim.B(~isnan(sim.B)));
if ( sumB > 6.7e5 && sumB < 6.8e5 )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

