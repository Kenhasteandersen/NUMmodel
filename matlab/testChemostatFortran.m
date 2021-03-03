function bSuccess = testChemostatFortran

sim = baserunChemostat(1.0, true);

sumB = sum(sim.B(:));
if ( sumB > 1430000 || sumB < 1400000 )
    bSuccess = false;
else
    bSuccess = true;
end

