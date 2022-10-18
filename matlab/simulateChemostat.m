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
    options.bUnicellularloss logical = true;
end
%
% Get the chemostat parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersChemostat(p, 'constantValues', [0.5 L]);
end
%
% Light
%
if ~isnan(p.seasonalOptions.lat_lon) | p.seasonalOptions.seasonalAmplitude~=0
    L = p.L;
end
sim.L = L;
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
uDeep(1) = p.uDeep; %Nutrients from the layer below the chemostat layer
%
% Check if there is POM:
%
ixGroupPOM = find(p.typeGroups==100);
if ~isempty(ixGroupPOM)
    ixPOM = p.ixStart(ixGroupPOM):p.ixEnd(ixGroupPOM);
else
    ixPOM = [];
end
p.velocity = 0*p.m;
p.velocity = calllib(loadNUMmodelLibrary(), 'f_getsinking', p.velocity);
%
% Simulate:
%
sLibname = loadNUMmodelLibrary();
[t,u] = ode45(@fDeriv, [0 p.tEnd], p.u0);
%
% Enforce minimum concentration
%
for k = 1:size(u,1)
    u(k,u(k,:)<p.umin) = p.umin(u(k,:)<p.umin);
end
%
% Enforce minimum concentration
%
for k = 1:size(u,1)
    u(k,u(k,:)<p.umin) = p.umin(u(k,:)<p.umin);
end

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
sim.rates = getRates(p, u(end,:), mean(sim.L), T);
for iGroup = 1:p.nGroups
    sim.Bgroup(:,iGroup) = sum( u(:, p.ixStart(iGroup):p.ixEnd(iGroup)),2);
end
sim.T = T;
sim.bUnicellularloss = options.bUnicellularloss;
%Bpnm = calcPicoNanoMicro(sim.B(end,:), sim.p.pGeneralists);
%sim.Bpico = Bpnm(1);
%sim.Bnano = Bpnm(2);
%sim.Bmicro = Bpnm(3);

[sim.Nbalance, sim.Cbalance,sim.Sibalance] = getBalance(sim.u(end,:), mean(sim.L), sim.T); % in units per day

    %
    % Function to assemble derivative for chemostat:
    %
    function dudt = fDeriv(t,u)
        dudt = 0*u';
        if (isnan(p.seasonalOptions.lat_lon) & p.seasonalOptions.seasonalAmplitude==0)
            [u, dudt] = calllib(sLibname, 'f_calcderivatives', ...
                u, L, T, 0.0, dudt);
            %
            % Chemostat dynamics for nutrients and unicellulars:
            %
            dudt(ix) = dudt(ix) + p.d*(uDeep(ix)-u(ix)');
        else % Incorporate the time dependency if necessary
            t_int = floor(mod(t,365))+1;
            if t_int>365
                t_int = 365;
            end
            [u, dudt] = calllib(sLibname, 'f_calcderivatives', ...
                u, L(t_int), T, 0.0, dudt);
            %
            % Chemostat dynamics for nutrients and unicellulars:
            %
            dudt(ix) = dudt(ix) + p.d(t_int)*(uDeep(ix)-u(ix)');
        end
        %
        % Sinking of POM:
        %
        dudt(ixPOM) = dudt(ixPOM) - p.velocity(ixPOM).*u(ixPOM)';
        dudt = dudt';
    end
end