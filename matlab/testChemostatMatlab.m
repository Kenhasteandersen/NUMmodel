function bSuccess = testChemostatMatlab

sim = baserunChemostat(1.0, false);

sumB = sum(sim.B(:));
if ( sumB > 930000 || sumB < 920000 )
    bSuccess = false;
else
    bSuccess = true;
end

