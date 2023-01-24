function bSuccess = testChemostatSeasonal(value)

p = setupGeneric([]);

p = parametersChemostat(p, 'lat_lon', [60 -10]);
sim = simulateChemostat(p, 'bUnicellularloss', false);

p_bis = parametersChemostat(p, 'seasonalAmplitude', 1);
sim_bis = simulateChemostat(p_bis, 'bUnicellularloss', false);

sumB = sum(sim.B(:));
sumB_bis = sum(sim_bis.B(:));
if ( sumB > 3963924  && sumB < 3963925 ) && ( sumB_bis > 1423361 && sumB_bis < 1423500 )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f; sum(B_bis) = %f\n", [sumB, sumB_bis]);
end