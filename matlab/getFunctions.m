function [ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro] = getFunctions(u, L)
%
% First make a call to calc a derivative:
%
u = double(u);
dudt = 0*u';
[u, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
            length(u), u, L, 0.0, dudt);
%
% Then extract the functions:
%
ProdGross = 0;
ProdNet = 0;
ProdHTL = 0;
eHTL = 0;
Bpico = 0;
Bnano = 0;
Bmicro = 0;

[ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro]...
    = calllib(loadNUMmodelLibrary(), 'f_getfunctions', ...
    ProdGross, ProdNet,ProdHTL,eHTL,Bpico,Bnano,Bmicro);
