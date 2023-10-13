function Bphyto_frac = calcBphytoFrac(jLreal,jF,jDOC,ix) %check units!
        Bphyto_frac = jLreal(ix)./(jLreal(ix)+jF(ix)+jDOC(ix));
end