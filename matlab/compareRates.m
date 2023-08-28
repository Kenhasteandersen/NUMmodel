%
% Calculate and compare jLreal, jFreal & jDOC from the 2 different simulations
%
% If tDay<0, average over the last month of the simulation
%
% In:
%   sim1 & sim2 - simulations to compare
%   tDay - time (-1 by default)
%

function Rates=compareRates(sim1,sim2,tDay)

arguments
    sim1 struct;
    sim2 struct;
    tDay double = -1;
end

% calculate the rates
rates(1) = calcRates(sim1,tDay);
rates(2) = calcRates(sim2,tDay);

% jLreal
Rates.jLreal = compareFields(rates(1).jLreal,rates(2).jLreal);

% jFreal
Rates.jFreal = compareFields(rates(1).jFreal,rates(2).jFreal);

% jDOC
Rates.jDOC = compareFields(rates(1).jDOC,rates(2).jDOC);
