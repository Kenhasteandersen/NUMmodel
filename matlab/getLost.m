%
% Get the losses of organic carbon, N, and Si
%
% In:
%  u - state variable vector (nutrients and biomasses of all groups).
%  L - light
%  T - temperature
%
% Out
%  Losses of organic carbon, N, and Si in units of XX/day
%
function [Clost, Nlost, SiLost] = getLost(u, L, T)
arguments
    u double;
    L double;
    T double;
end
%
% Make a call to calc derivatives:
%
u = double(u);
dudt = 0*u';
[u, dudt]=calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
            u, L, T, 0.0, dudt);
%
% Then extract losses:
%
Clost = 0;
Nlost = 0;
SiLost= 0;

[~, Clost, Nlost, SiLost] = calllib(loadNUMmodelLibrary(), 'f_getlost', ...
    u, Clost, Nlost, SiLost);