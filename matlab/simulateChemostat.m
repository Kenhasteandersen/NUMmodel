
%
% Simulate the chemostat.  To run a seasonally varying simulation, set the
% "lat_lon" option in "parametersChemostat".
%
% In:
%  p - parameter object obtained by calling a "setup" function followed
%      by a call to parametersChemostat; see e.g. the default value below.
%  L - Light
%  T - Temperature
% Options:
%  bUnicellularloss - determines whether unicellular groups are subject to
%      mixing losses
%  bVerbose - displays the nutrients balances
%
% Out:
%  sim - simulation object. The N budget of the layer is integrated along
%        with the state (units mugN/l, cumulative from t=0):
%          Nprod    - N mixed in from the deep layer
%          Nloss    - N lost from the layer: mixing of unicellulars, sinking
%                     of POM, and (without a POM group) HTL and POM export
%          NlossHTL - the HTL/POM-export part of Nloss (zero with POM)
%
function sim = simulateChemostat(p, L, T, options)

arguments
    p struct = parametersChemostat(setupGeneralistsOnly);
    L double = 100;
    T double = 10;
    options.bUnicellularloss logical = true;
    options.bVerbose logical = false;
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
    ix = 1:p.ixEnd( max(find(p.typeGroups<10)) ); % Nitrogen and unicellulars are lost to the deep layer
else
    ix = 1:(p.idxB-1); % Nutrients
end
ixUniB = ix(ix>=p.idxB); % The biomass part of what is mixed
uDeep = p.uDeep;
uDeep(p.idxB:length(p.u0)) = 0;
bSeasonal = ~isnan(p.seasonalOptions.lat_lon) | p.seasonalOptions.seasonalAmplitude~=0;
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

rhoCN = search_namelist('../input/input.yaml','general','rhoCN');
% The state is augmented with the three cumulative N-budget terms:
[t,y] = ode23s(@fDeriv, [0 p.tEnd], [p.u0 0 0 0]);
u = y(:,1:p.n);
sim.Nprod = y(:,p.n+1);
sim.Nloss = y(:,p.n+2);
sim.NlossHTL = y(:,p.n+3);

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

%
% Get the balance of the derivative:
%
[sim.Cbalance,sim.Nbalance,sim.Sibalance] = getBalance(sim.u(end,:), mean(sim.L), sim.T); % in units per day

%
% Display rates
%
if options.bVerbose
    %Rate
    Crate=sim.Cbalance/(sim.DOC(end)+sum(sim.B(end,:)))*100;
    Nrate=sim.Nbalance/(sim.N(end)+sum(sim.B(end,:))/rhoCN)*100;
    fprintf("----------------------------------------------\n")
    fprintf("Rate of gain of C: %8.3f %% per day \n", Crate);
    fprintf("Rate of gain of N in the last derivative calculation: %8.3f  %% per day \n", Nrate);
    %Presence of Diatoms 
    ixDiatoms = find(p.typeGroups==3);
    if ~isempty(ixDiatoms)
        ixDiatoms = (p.ixStart(ixDiatoms):p.ixEnd(ixDiatoms))-p.idxSi;
        rhoCSi = search_namelist('../input/input.yaml','diatoms','rhoCSi');
        Sirate=sim.Sibalance/(sim.Si(end)+sum(sim.B(end,ixDiatoms))/rhoCSi)*100;
        fprintf("Rate of gain of Si: %8.3f %% per day \n", Sirate);
    end
    %N budget over the simulation:
    fprintf("N mixed in from the deep: %8.3f mugN/l; N lost from the layer: %8.3f mugN/l\n", ...
        sim.Nprod(end), sim.Nloss(end));
    fprintf("----------------------------------------------\n")
end


    % -------------------------------------------------------------------------
    % Function to assemble derivative for chemostat:
    %
    function dydt = fDeriv(t,y)

        u = y(1:p.n)';
        if bSeasonal
            t_int = min(floor(mod(t,365))+1, 365);
            Lnow = L(t_int);
            dnow = p.d(t_int);
        else
            Lnow = L;
            dnow = p.d;
        end
        dudt = 0*u;
        [u, dudt] = calllib(sLibname, 'f_calcderivatives', ...
            u, Lnow, T, 0.0, dudt);
        %
        % Chemostat dynamics for nutrients and unicellulars:
        %
        dudt(ix) = dudt(ix) + dnow*(uDeep(ix)-u(ix));
        %
        % Sinking of POM:
        %
        sinkPOM = p.velocity(ixPOM).*u(ixPOM)/p.widthProductiveLayer;
        dudt(ixPOM) = dudt(ixPOM) - sinkPOM;
        %
        % N budget of the layer. The export via HTL and POM comes from the
        % library (zero when a POM group is present):
        %
        Nlost = 0;
        [~,~,Nlost,~] = calllib(sLibname, 'f_getlost', u, 0, Nlost, 0);
        Nprod = dnow*(uDeep(p.idxN)-u(p.idxN));
        Nloss = dnow*sum(u(ixUniB))/rhoCN + sum(sinkPOM)/rhoCN + Nlost;

        dydt = [dudt, Nprod, Nloss, Nlost]';
    end
end