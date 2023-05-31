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
%  bCalculateNloss - determines whether to calculate the variation of Nitrogen per day in the Chemostat layer (takes time)
%                    Adds to the output : NDeepTot - N gain from the deep throughout the simulation in µg/L
%                                         deltaNdt - Variation of Nitrogen per day in the Chemostat layer in µg/L/day
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
    options.bCalculateNloss logical = false;
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
if options.bCalculateNloss
    u0=[p.u0 0];
else 
    u0=p.u0;
end

sim.rhoCN=5.68;
[t,u] = ode23s(@fDeriv, [0 p.tEnd], u0); 

if options.bCalculateNloss
    sim.NDeepTot=u(:,end); % N gain from the deep throughout the simulation
    u=u(:,1:p.n);
end
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

%
% Calculate N variations
%
if options.bCalculateNloss
    sim.deltaNdt = (sim.N(end)+sum(sim.B(end,:))/sim.rhoCN - ...
        (p.u0(p.idxN)+sum(p.u0(p.idxB:end))/sim.rhoCN + ...
        sim.NDeepTot(end)))/ t(end); % Variation of Nitrogen per day in the Chemostat layer
end

[sim.Cbalance,sim.Nbalance,sim.Sibalance] = getBalance(sim.u(end,:), mean(sim.L), sim.T); % in units per day



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
    
            %
            % Calculate N gain from the deep
            %
            if options.bCalculateNloss
                
                % Extract the losses
                Clost=0;
                Nlost = 0;
                SiLost=0;

                [~,~, Nlost, ~] = calllib(sLibname, 'f_getlost', ...
                    u, Clost, Nlost, SiLost);
           
                dudt(end+1) = (uDeep(1)-u(1))*p.d-Nlost; 
            
                if options.bUnicellularloss 
                    dudt(end)=dudt(end)-p.d*sum(u(p.idxB:end))/sim.rhoCN; %takes B's losses to the deep into account
                end
            end

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
        dudt(ixPOM) = dudt(ixPOM) - p.velocity(ixPOM).*u(ixPOM)'/p.widthProductiveLayer;
        dudt = dudt';
    end
end