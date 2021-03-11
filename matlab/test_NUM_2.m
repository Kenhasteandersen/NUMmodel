clearvars

mAdult=10;

sim=baserunChemostat_csp(mAdult, true);

% Compare calculated rates:
unloadlibrary(loadNUMmodelLibrary());
loadNUMmodelLibrary();

%would be nice to get the initial conditions we used in fortran without checking them 
%(or simply state them from matlab))
u=zeros(size(sim.rates.dudt));
u(1) = 150.; %N
u(2) = 1.; %DOM
u(3:end) = 1e-2; %Plankton
dudt = 0*u;
L = 150;

% calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly');
calllib(loadNUMmodelLibrary(), 'f_setupgeneric_csp', int32(length(mAdult)), mAdult);
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', length(u), u, L, 0.0, dudt);
ratesF = calcRates(u,L);

p = parameters([]);
% rates = calcDerivatives(p, u, L);

plot_diagnostics(p,sim,ratesF)

% disp('jN:')
% [ratesF.jN'; rates.JN./p.m]
% 
% disp('jL:')
% [ratesF.jL'; rates.JL./p.m]
% 
% disp('jF:')
% [ratesF.jF'; rates.JF./p.m]
% 
% disp('jTot:')
% [ratesF.jTot'; rates.Jtot./p.m]
% 
% disp('mortHTL:')
% [ratesF.mortHTL'; p.mortHTLm]
% 
% disp('mortpred:')
% [ratesF.mortpred'; rates.mortpred]
% 
% disp('g:')
% [ratesF.g'; rates.JF./p.m]
% 
% 
% disp('dudt:')
% [dudt; rates.dudt]


