function bSuccess = testChemostatMatlab

sim = baserunChemostat(1.0, false);

sumB = sum(sim.B(:));
if ( sumB > 7e6 || sumB < 6.7e6 )
    bSuccess = false;
else
    bSuccess = true;
end

