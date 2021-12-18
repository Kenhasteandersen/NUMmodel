function bSuccess = testChemostat

sim = baserunChemostat(1.0);

sumB = sum(sim.B(:));
if ( sumB > 2.2e7 || sumB < 2e7 )
    bSuccess = false;
else
    bSuccess = true;
end

