
%clearvars

%baserunChemostat([], true)

%%
% Compare calculated rates:

unloadlibrary(loadNUMmodelLibrary());
loadNUMmodelLibrary();
n = 10;

u(1) = 150.;
u(2) = 1.;
u(3:n+2) = 1:n;
%u(13:22) = 1.;
%u(3:22) = 1:20;
dudt = 0*u;
L = 150;

%calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(1), 0.1);
%calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly', int32(n));
calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly_csp');
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', length(u), u, L, 0.0, dudt);
ratesF = calcRates(u,L);

p = parameters([]);
rates = calcDerivatives(p, u, L);

disp('jN:')
[ratesF.jN'; rates.JN./p.m]

disp('jL:')
[ratesF.jL'; rates.JL./p.m]

disp('jF:')
[ratesF.jF'; rates.JF./p.m]

disp('jTot:')
[ratesF.jTot'; rates.Jtot./p.m]

disp('mortHTL:')
[ratesF.mortHTL'; p.mortHTLm]

disp('mortpred:')
[ratesF.mortpred'; rates.mortpred]

disp('g:')
[ratesF.g'; rates.JF./p.m]


disp('dudt:')
[dudt; rates.dudt]
