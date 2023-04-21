%
% Get the ecosystem functions.
%
% In:
%  u - state variable vector (nutrients and biomasses of all groups).
%  L - light
%  T - temperature
%
% Out
%  ProdGross, ProdNet,ProdHTL - Gross, net, and HTL productions in units of
%                               gC per m^3.
%  eHTL - not implemented
%  Bpico, Bnano, Bmicro - biomasses in pico, nano, and micro size groups
%                         (gC/m3).
%
function [ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro] = ...
    getFunctions(u, L, T, sLibName)
arguments
    u double;
    L double;
    T double;
    sLibName = loadNUMmodelLibrary();
end

%
% First make a call to calc a derivative:
%
u = double(u);
dudt = 0*u';
[u, dudt] = calllib(sLibName, 'f_calcderivatives', ...
            u, L, T, 0.0, dudt);
%
% Then extract the functions:
%
ProdGross = 0;
ProdNet = 0;
ProdHTL = 0;
ProdBact = 0;
eHTL = 0;
Bpico = 0;
Bnano = 0;
Bmicro = 0;

[u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro]...
    = calllib(sLibName, 'f_getfunctions', ...
    u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro);
