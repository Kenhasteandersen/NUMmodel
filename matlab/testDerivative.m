p = parameters(100);

rates = calcDerivatives(p,p.u0,100);

plotRates(p,rates)
