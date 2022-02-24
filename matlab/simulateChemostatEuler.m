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
function sim = simulateChemostatEuler(p, L, T, bLosses)

arguments
    p struct = parametersChemostat(setupGeneralistsOnly);
    L double = 100;
    T double = 10;
    bLosses logical = false;
end
%
% Get the chemostat parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersChemostat(p);
end

%
% Simulate:
%
u = p.u0;
dudt = 0*u;

u = calllib(loadNUMmodelLibrary(), 'f_simulatechemostateuler', u, ...
    L, T, ...
    int32(p.idxB-1), ...
    p.u0(1:(p.idxB-1)), p.d, p.tEnd, 0.01, bLosses);
%
% Assemble result:
%
sim.t = p.tEnd;
sim.N = u(p.idxN);
sim.DOC = u(p.idxDOC);
if isfield(p, 'idxSi')
    sim.Si = u(p.idxSi);
end
sim.B = u(3:end);
sim.p = p;
sim.rates = getRates(sim.p, u(end,:), L, T);
for iGroup = 1:p.nGroups
    sim.Bgroup(:,iGroup) = sum( u(p.ixStart(iGroup):p.ixEnd(iGroup)));
end
sim.L = L;
sim.T = T;

[sim.Nbalance, sim.Cbalance] = getBalance(u, sim.L, sim.T); % in units per day

end