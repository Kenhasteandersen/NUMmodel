function bSuccess = testChemostatFortran

sim = baserunChemostat(1.0, true);

sumB = sum(sim.B(:));
if ( sumB > 2.2e7 || sumB < 2e7 )
    bSuccess = false;
else
    bSuccess = true;
end

