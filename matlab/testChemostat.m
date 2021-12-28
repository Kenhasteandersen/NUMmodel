function bSuccess = testChemostat

sim = baserunChemostat(1.0);

sumB = sum(sim.B(:));
if ( sumB > 3.9e7 && sumB < 4e7 )
    bSuccess = true;
else
    bSuccess = false;
end

