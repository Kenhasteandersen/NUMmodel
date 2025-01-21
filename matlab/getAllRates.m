%
% Get all the rates in a simulation.
% Only works with chemostat.
%
% In:
%  sim - the simulation structure
%
% Out:
%  rates - a matrix with the rates at all time steps
%
function rates = getAllRates(sim)

p = sim.p;
if ~strcmp(p.nameModel,'chemostat')
    error('Only works on chemostat');
end

for i = 1:length(sim.t)
    % Get light and tempertaure:
    if length(sim.L)>1
        L = interp1(1:length(sim.L), sim.L, sim.t(i) );
    else
        L = sim.L;
    end
    T = sim.T;
    % Assemble the state vector:
    if p.nNutrients==3
        u = [sim.N(i), sim.DOC(i),sim.Si(i), sim.B(i,:)];
    else
        u = [sim.N(i), sim.DOC(i), sim.B(i,:)];
    end

    rate = getRates(p, u, L, T);

    rates.jN(i,:) = rate.jN;
    rates.jDOC(i,:) = rate.jDOC;
    rates.jL(i,:) = rate.jL;
    rates.jSi(i,:) = rate.jSi;
    rates.jF(i,:) = rate.jF;
    rates.jFreal(i,:) = rate.jFreal;
    rates.f(i,:) = rate.f;
    rates.jTot(i,:) = rate.jTot;
    rates.jMax(i,:) = rate.jMax;
    rates.jFmaxx(i,:) = rate.jFmaxx;
    rates.jR(i,:) = rate.jR;
    rates.jRespTot(i,:) = rate.jRespTot;
    rates.jLossPassive(i,:) = rate.jLossPassive;
    rates.jNloss(i,:) = rate.jNloss;
    rates.jLreal(i,:) = rate.jLreal;
    rates.jPOM(i,:) = rate.jPOM;
    rates.mortpred(i,:) = rate.mortpred;
    rates.mortHTL(i,:) = rate.mortHTL;
    rates.mort2(i,:) = rate.mort2;
    rates.mort(i,:) = rate.mort;
end

