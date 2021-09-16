%
% Returns all rates from the unicellulars for a given state (u) and light
% level.
%
function rates = getRates(p, u, L, T)
arguments
    p struct;
    u double;
    L double;
    T double = 10;
end
%
% First make a call to calc a derivative:
%
dudt = 0*u';
calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
            length(u), u, L, T, 0.0, dudt);
%
% Then extract the rates:
%
zero = zeros(length(u)-p.idxB+1,1);
jN = zero;
jDOC = zero; 
jL = zero; 
jSi = zero;
jF = zero;
jFreal = zero;
jTot = zero;
jMax = zero;
jFmaxx = zero;
jR = zero;
jLossPassive = zero;
jNloss = zero;
jLreal = zero;
mortpred = zero;
mortHTL = zero;
mort2 = zero;
mort = zero;

[rates.jN, rates.jDOC, rates.jL, rates.jSi, rates.jF, rates.jFreal, rates.jTot, rates.jMax, rates.jFmaxx, ...
    rates.jR, rates.jLossPassive, rates.jNloss, rates.jLreal, rates.mortpred, rates.mortHTL, ...
    rates.mort2, rates.mort] = ...
    calllib(loadNUMmodelLibrary(), 'f_getrates', ...
    jN, jDOC, jL, jSi, jF, jFreal,...
    jTot, jMax, jFmaxx, jR, jLossPassive, ...
    jNloss,jLreal, ...
    mortpred, mortHTL, mort2, mort);


