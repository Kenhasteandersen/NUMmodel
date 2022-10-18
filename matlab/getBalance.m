%
% Get the nitrogen and carbon balance. Only works with setupGeneralistsOnly
%
% In:
%  u - state variable vector (nutrients and biomasses of all groups).
%  L - light
%  T - temperature
%
% Out
%  N balance, C balance in units 1/day

%
function [Nbalance, Cbalance,Sibalance] = getBalance(u, L, T)
arguments
    u double;
    L double;
    T double;
end
u = double(u);
dudt = 0*u';
[u, dudt]=calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
            u, L, T, 0.0, dudt);
%
% Then extract balance:
%
Nbalance = 0;
Cbalance = 0;
Sibalance=0;

[~, ~, Nbalance, Cbalance,Sibalance] = calllib(loadNUMmodelLibrary(), 'f_getbalance', ...
    u, dudt, Nbalance, Cbalance, Sibalance);

