%
% Returns rates calculated from Fortran library
%
function rates = calcRates(u,L)
loadNUMmodelLibrary();

jN = 0*u';
jL = jN;
jF = jN;

[u, jN,jL,jF] = calllib(loadNUMmodelLibrary(), 'f_calcrates', ...
            length(u), u, L, jN,jL,jF);

rates.jN = jN;
rates.jL = jL;
rates.jF = jF;