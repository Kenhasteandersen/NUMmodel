%
% Get the carbon, nitrogen and silicate balances in the calculation of the
% derivatives.
%
% In:
%  u - state variable vector (nutrients and biomasses of all groups).
%  L - light
%  T - temperature
%
% Out
%  N balance, C balance, and Si balance in units 1/day

%
function [Cbalance, Nbalance, Sibalance, dudt] = getBalance(u, L, T, bPrintSummary)
arguments
    u double;
    L double;
    T double;
    bPrintSummary = false;
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
    u, dudt, Cbalance, Nbalance, Sibalance);

if bPrintSummary
    fprintf("----------------------------------------------\n")
    fprintf(" Conservation of carbon    %6.2e 1/day\n", Cbalance); 
    fprintf(" Conservation of nutrients %6.2e 1/day\n", Nbalance); 
    fprintf(" Conservation of silicate  %6.2e 1/day\n", Sibalance); 
    fprintf("----------------------------------------------\n")
end

