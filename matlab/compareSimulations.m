% 
% Compare 2 simulations structures with same parameters
%
% In:
%   sim1 & sim2 - simulation structures to compare
%   Optional:
%   options.rates - compare jLreal, jFreal & jDOC from the 2 simulations 
%
% Out:
%   sim - difference between sim1 et sim2 in %
%

function sim=compareSimulations(sim1,sim2,options)

arguments
    sim1 struct;
    sim2 struct;
    options.rates logical = false;
end

sim=sim1;

sim.N = compareFields(sim1.N,sim2.N);
sim.DOC = compareFields(sim1.DOC,sim2.DOC);
if isfield(sim,'Si')
    sim.Si = compareFields(sim1.Si,sim2.Si);
end
sim.B = compareFields(sim1.B,sim2.B);
sim.L = compareFields(sim1.L,sim2.L);
sim.T = compareFields(sim1.T,sim2.T);

if options.rates
    sim.optionRates = options.rates;
    sim.rates=compareRates(sim1,sim2,-1);
end













