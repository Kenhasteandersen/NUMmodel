function bSuccess = testChemostat(value)

sim = baserunChemostat([]);

sumB = sum(sim.B(:));
if ( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

