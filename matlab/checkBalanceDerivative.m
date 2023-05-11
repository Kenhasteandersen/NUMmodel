%
% Checks for N balance in the calculation of the derivative from the
% Fortran NUMmodel library.  It does so by calculating the derivative for
% the last time step in a simulation.
%
% In:
%   sim - a simulation object to check
%
% Out:
%   The amount of N unaccounted for
%
function Nbalance = checkBalanceDerivative(sim)
%
% Constants:
%
fracHTL_to_N = search_namelist('../input/input.nlm','general','fracHTL_to_N');
rhoCN = search_namelist('../input/input.nlm','general','rhoCN');
%
% Extract u, dudt, and rates from the last time step:
%
p = sim.p;
if ~isfield(sim.p,'idxSi')
    u = [sim.N(end), sim.DOC(end), sim.B(end,:)];
else
    u = [sim.N(end), sim.DOC(end), sim.Si(end), sim.B(end,:)];
end

dudt = 0*u;
[~, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
                u, sim.L(end), sim.T(end), 0.0, dudt);
rates = getRates(sim.p,u,sim.L(end), sim.T(end));

B = u(p.idxB:end);
if ~sum(ismember(p.typeGroups,100))
    % Losses from HTL:
    lossHTL = sum((1-fracHTL_to_N)*rates.mortHTL.*B')/rhoCN;
    lossPOM = sum(rates.jPOM.*B') / rhoCN;
else
    lossHTL = 0;
    ixPOM = p.ixStart(p.ixPOM):p.ixEnd(p.ixPOM);
    lossPOM = p.velocity(ixPOM).*u(ixPOM)/rhoCN;
end

Nbalance = dudt(p.idxN) + sum(dudt(p.idxB:end))/rhoCN + lossHTL + lossPOM;