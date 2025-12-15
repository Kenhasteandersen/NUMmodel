function bSuccess = testChemostat(value)

p = setupNUMmodel();       % Sets up the model
p = parametersChemostat(p);% Sets up the chemostat environment
sim = simulateChemostat();

sumB = sum(sim.B(:));
if ( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

