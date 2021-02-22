%
% Returns rates calculated from Fortran library
%
function rates = calcRates(u,L)
loadNUMmodelLibrary();

jN = 0*u';
jL = jN;
jF = jN;
jTot = jN;
mortHTL = jN;
mortpred = jN;
g = jN;

[u, jN,jL,jF,jTot,mortHTL,mortpred,g] = calllib(loadNUMmodelLibrary(), 'f_calcrates', ...
            length(u), u, L, jN,jL,jF,jTot,mortHTL,mortpred,g);

rates.jN = jN;
rates.jL = jL;
rates.jF = jF;
rates.jTot = jTot;
rates.mortHTL = mortHTL;
rates.mortpred = mortpred;
rates.g = g;