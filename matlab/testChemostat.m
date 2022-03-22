function bSuccess = testChemostat

sim = baserunChemostat(1.0);

sumB = sum(sim.B(:));
if ( sumB > 1.27e7 && sumB < 1.28e7 )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

