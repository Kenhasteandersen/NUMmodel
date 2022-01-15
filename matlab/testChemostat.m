function bSuccess = testChemostat

sim = baserunChemostat(1.0);

sumB = sum(sim.B(:));
if ( sumB > 4.8e7 && sumB < 4.9e7 )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

