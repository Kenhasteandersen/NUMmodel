
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
%  bCalculateNgain - determines whether to calculate the variation of Nitrogen per day in the Chemostat layer
%                    Adds to the output : Ngain - N gained from the deep throughout the simulation minus the losses in µg/L
%                                         deltaNdt - Variation of N per day in the Chemostat layer in µg/L/day
%  bVerbose - displays the nutrients balances
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
    options.bCalculateNgain logical = false;
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
uDeep = p.uDeep;
uDeep(p.idxB:length(p.u0)) = 0;
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

% Initial conditions
if options.bCalculateNgain
    u0=[p.u0 0];
else 
    u0=p.u0;
end

rhoCN = search_namelist('../input/input.h','general','rhoCN');
[t,u] = ode23s(@fDeriv, [0 p.tEnd], u0); 

if options.bCalculateNgain
    sim.Ngain=u(:,end); % N gained from the deep minus the losses
    u=u(:,1:p.n);
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

%
% Calculate N variations
%
if options.bCalculateNgain
    sim.deltaNdt = (sim.N(end)+sum(sim.B(end,:))/rhoCN - ...
        (p.u0(p.idxN)+sum(p.u0(p.idxB:end))/rhoCN + ...
        sim.Ngain(end)))/ t(end); % Variation of Nitrogen per day in the Chemostat layer
end
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
        rhoCSi = search_namelist('../input/input.h','diatoms','rhoCSi');
        Sirate=sim.Sibalance/(sim.Si(end)+sum(sim.B(end,ixDiatoms))/rhoCSi)*100;
        fprintf("Rate of gain of Si: %8.3f %% per day \n", Sirate);
    end
    %N losses
    if options.bCalculateNgain
        changes=sim.deltaNdt/(sim.N(end)+sum(sim.B(end,:))/5.68)*100;
        fprintf("Total gain of N throughout the simulation: %8.3f µg/L \n", sim.Ngain(end));
        fprintf("Average rate of N change over the entire simulation: %8.3f %% per day \n", changes)
    end
    fprintf("----------------------------------------------\n")
end


    % -------------------------------------------------------------------------
    % Function to assemble derivative for chemostat:
    %
    function dudt = fDeriv(t,u)

        u = u(1:p.n);
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
        % Calculate N gain from the deep
        %
        if options.bCalculateNgain
            
            % Extract the losses
            Clost=0;
            Nlost = 0;
            SiLost=0;

            [~,~, Nlost, ~] = calllib(sLibname, 'f_getlost', ...
                u, Clost, Nlost, SiLost);
       
            dudt(end+1) = (uDeep(1)-u(1))*p.d-Nlost;
        
            if options.bUnicellularloss 
                dudt(end)=dudt(end)-p.d*sum(u(p.idxB:end))/rhoCN; %takes B's losses to the deep into account
            end
        end
        %
        % Sinking of POM:
        %
        dudt(ixPOM) = dudt(ixPOM) - p.velocity(ixPOM).*u(ixPOM)'/p.widthProductiveLayer;
        dudt = dudt';
    end
end