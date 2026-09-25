%
% Simulate the chemostat with Euler integration in Fortran library
% In:
%  p - parameter object (including chemostat parameters from
%      parametersChemostat). Not used if running from the fortran library. 
%  L - Light
%  T - Temperature
%  bLosses - whether to have losses to the deep
%
% Out:
%  sim - simulation object
%
function sim = simulateChemostatEuler(p, L, T, options)

arguments
    p struct = parametersChemostat(setupGeneralistsOnly);
    L double = 100;
    T double = 10;
    options.bUnicellularloss logical = true;
end
%
% Get the chemostat parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersChemostat(p);
end
%
% Simulate. The mixing with the deep layer is handled inside the library:
%
u = p.u0;

u = calllib(loadNUMmodelLibrary(), 'f_simulatechemostateuler', u, ...
    L, T, ...
    int32(p.idxB-1), ...
    p.uDeep(1:(p.idxB-1)), p.d, p.widthProductiveLayer, ...
    p.tEnd, 0.01, options.bUnicellularloss);
%
% Functions of the solution. Note that f_simulateeulerfunctions must not be
% used here: it integrates u for a further tEnd days without the chemostat
% dynamics before evaluating the functions.
%
[sim.ProdGross, sim.ProdNet, sim.ProdHTL, sim.ProdBact, sim.eHTL, ...
    sim.Bpico, sim.Bnano, sim.Bmicro, sim.mHTL] = getFunctions(u, L, T);


%
% Assemble result:
%
sim.t = p.tEnd;
sim.u = u;
sim.N = u(p.idxN);
sim.DOC = u(p.idxDOC);
if isfield(p, 'idxSi')
    sim.Si = u(p.idxSi);
end
sim.B = u(p.idxB:end);
sim.p = p;
sim.rates = getRates(sim.p, u(end,:), L, T);
for iGroup = 1:p.nGroups
    sim.Bgroup(:,iGroup) = sum( u(p.ixStart(iGroup):p.ixEnd(iGroup)));
end
sim.L = L;
sim.T = T;
sim.bUnicellularloss = options.bUnicellularloss;

[sim.Cbalance, sim.Nbalance, sim.Sibalance] = getBalance(u, sim.L, sim.T); % in units per day

end