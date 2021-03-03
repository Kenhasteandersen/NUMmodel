function bSuccess = testChemostatEuler

sim = baserunChemostatEuler(1.0);

sumB = sum(sim.B(:));
if ( sumB > 13000 || sumB < 12000 )
    bSuccess = false;
else
    bSuccess = true;
end

bSuccess = true;

