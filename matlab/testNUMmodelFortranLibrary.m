% unloadlibrary(loadNUMmodelLibrary());
% loadNUMmodelLibrary();
% 
% u(1) = 150.;
% u(2) = 1.;
% u(3:12) = 1:10;
% u(13:22) = 1.;
% u(3:22) = 1:20;
% dudt = 0*u;
% L = 100;
% calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(1), 0.1);
% [u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', length(u), u, L, 0.0, dudt);
% 
% 
% p = parameters(.1);
% rates = calcDerivatives(p, u, L);
% 
% format longe
% %rates.F./p.m
% 
% rates.dudt
% dudt


%%

clearvars

 baserunChemostat([0.1 1 10 100 1000])