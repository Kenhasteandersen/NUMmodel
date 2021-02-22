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
g = jN;

[u, jN,jL,jF,jTot,mortHTL,g] = calllib(loadNUMmodelLibrary(), 'f_calcrates', ...
            length(u), u, L, jN,jL,jF,jTot,mortHTL,g);

rates.jN = jN;
rates.jL = jL;
rates.jF = jF;
rates.jTot = jTot;
rates.mortHTL = mortHTL;
rates.g = g;