function bSuccess = testChemostatSeasonal(value)

p = setupNUMmodel();

p = parametersChemostat(p, 'lat_lon', [60 -10]);
sim = simulateChemostat(p, 'bUnicellularloss', false);

p_bis = parametersChemostat(p, 'seasonalAmplitude', 1);
sim_bis = simulateChemostat(p_bis, 'bUnicellularloss', false);

sumB = sum(sim.B(:));
sumB_bis = sum(sim_bis.B(:));
if ( sumB > 1479637  && sumB < 1479638 ) && ( sumB_bis > 24752 && sumB_bis < 24753 )
    bSuccess = true;
else
    bSuccess = false;
    fprintf(2,"sum(B) = %f; sum(B_bis) = %f\n", [sumB, sumB_bis]);
end