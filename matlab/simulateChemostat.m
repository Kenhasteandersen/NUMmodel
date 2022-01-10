%
% Simulate the chemostat.
% In:
%  p - parameter object obtained by calling a "setup" function followed
%      by a call to parametersChemostat; see e.g. the default value below.
%  L - Light
%  T - Temperature
% Options:
%  bUnicellularloss - determines whether unicellular groups are subject to
%      mixing losses
%
% Out:
%  sim - simulation object
%
function sim = simulateChemostat(p, L, T, options)

arguments
    p struct = parametersChemostat(setupGeneralistsOnly);
    L double = 100;
    T double = 10;
    options.bUnicellularloss logical = false;
end
%
% Concentrations in the deep layer:
%
if options.bUnicellularloss
    ix = 1:p.ixEnd(1); % Nutrients and first field (unicellulars) are lost to the deep layer
else
    ix = 1:(p.idxB-1); % Nutrients
end
uDeep = p.u0;
uDeep(p.idxB:end) = 0;
%
% Simulate:
%
sLibname = loadNUMmodelLibrary();
[t,u] = ode23(@fDeriv, [0 p.tEnd], p.u0);
%
% Assemble result:
%
sim.u=u;
sim.t = t;
sim.N = u(:,p.idxN);
sim.DOC = u(:,p.idxDOC);
if isfield(p, 'idxSi')
    sim.Si = u(:,p.idxSi);
end
sim.B = u(:,p.idxB:end);
sim.p = p;
sim.rates = getRates(p, u(end,:), L);
for iGroup = 1:p.nGroups
    sim.Bgroup(:,iGroup) = sum( u(:, p.ixStart(iGroup):p.ixEnd(iGroup)),2);
end
sim.L = L;
sim.T = T;
%Bpnm = calcPicoNanoMicro(sim.B(end,:), sim.p.pGeneralists);
%sim.Bpico = Bpnm(1);
%sim.Bnano = Bpnm(2);
%sim.Bmicro = Bpnm(3);

[sim.Nbalance, sim.Cbalance] = getBalance(sim.u, sim.L, sim.T); % in units per day

    %
    % Function to assemble derivative for chemostat:
    %
    function dudt = fDeriv(t,u)
        dudt = 0*u';
        [u, dudt] = calllib(sLibname, 'f_calcderivatives', ...
            u, L, T, 0.0, dudt);
        %
        % Chemostat dynamics for nutrients and unicellulars:
        %
        dudt(ix) = dudt(ix) + p.d*(uDeep(ix)-u(ix)');
        dudt = dudt';
    end
end