function bSuccess = testChemostat(value)

p = setupNUMmodel();       % Sets up the model
p = parametersChemostat(p);% Sets up the chemostat environment
sim = simulateChemostat();

bSuccess = true;

sumB = sum(sim.B(:));
if ~( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

[~, dNdt_per_N] = checkConservation(sim, false);
fprintf("relative N balance = %+.2e 1/yr. ", dNdt_per_N);
if abs(dNdt_per_N) > 1e-4
    bSuccess = false;
    fprintf(2,"N not conserved!\n");
end
