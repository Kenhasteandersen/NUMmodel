function bSuccess = testChemostat

sim = baserunChemostat(1.0);

sumB = sum(sim.B(:));
if ( sumB > 4.e7 && sumB < 4.1e7 )
    bSuccess = true;
else
    bSuccess = false;
end

