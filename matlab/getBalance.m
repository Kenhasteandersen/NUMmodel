%
% Get the carbon, nitrogen and silicate balances.
%
% In:
%  u - state variable vector (nutrients and biomasses of all groups).
%  L - light
%  T - temperature
%
% Out
%  N balance, C balance, and Si balance in units 1/day

%
function [Cbalance, Nbalance, Sibalance] = getBalance(u, L, T)
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
Cbalance = 0;
Nbalance = 0;
Sibalance= 0;

[~, ~,Cbalance, Nbalance,Sibalance] = calllib(loadNUMmodelLibrary(), 'f_getbalance', ...
    u, dudt, Cbalance, Nbalance,Sibalance);
