function bSuccess = testWatercolumn(value)

p = setupGeneralistsPOM(5,1); % A fast setup with POM
p = parametersWatercolumn(p);

sim = simulateWatercolumn(p,60,-15);
checkConservation(sim);

sumB = sum(sim.B(~isnan(sim.B)));
if ( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end
