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
% Concentrations in the deep layer:
%
if options.bUnicellularloss
    ix = 1:p.ixEnd(1); % Nutrients and first field (unicellulars) are lost to the deep layer
else
    ix = 1:(p.idxB-1); % Nutrients
end
uDeep = p.uDeep;
uDeep(p.idxB:length(p.u0)) = 0;
%
% Simulate:
%
u = p.u0;
dudt = 0*u;

u = calllib(loadNUMmodelLibrary(), 'f_simulatechemostateuler', u, ...
    L, T, ...
    int32(p.idxB-1), ...
    p.u0(1:(p.idxB-1)), p.d, p.tEnd, 0.01, options.bUnicellularloss);

ProdGross = 0;
ProdNet = 0;
ProdHTL = 0;
ProdBact = 0;
eHTL = 0;
Bpico = 0;
Bnano = 0;
Bmicro = 0;
[u, ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro] = ...
    calllib(loadNUMmodelLibrary(), 'f_simulateeulerfunctions', u, ...
    L, T, p.tEnd, 0.01,  ...
    ProdGross, ProdNet,ProdHTL,ProdBact,eHTL,Bpico,Bnano,Bmicro);


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