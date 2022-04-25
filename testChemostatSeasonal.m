function bSuccess = testChemostatSeasonal

p = setupGeneric([]);

p = parametersChemostat(p, 'lat_lon', [60 -10]);
sim = simulateChemostat(p, 'bUnicellularloss', false);

p_bis = parametersChemostat(p, 'seasonalAmplitude', 1);
sim_bis = simulateChemostat(p_bis, 'bUnicellularloss', false);

sumB = sum(sim.B(:));
sumB_bis = sum(sim_bis.B(:));
if ( sumB > 2e5 && sumB < 3e5 ) && ( sumB_bis > 9e4 && sumB_bis < 10e4 )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end