%
% Calculates the nitrogen balance of the course of the simulation
%
% In:
%  sim: simulation structure
%  bVerbose: whether to print out the balance on the terminal (default true
%           ).
%
function [dNdt, dNdt_per_N] = checkConservation(sim, bVerbose)

arguments
    sim struct;
    bVerbose logical = true; % Whether to print out the solution
end
%
% Constants:
%
epsilonF = 0.8;
reminHTL = 0;
reminF = 0.1;
rhoCN = 5.68;

p = sim.p;
gains = 0;
losses = 0;
dt = gradient(sim.t);

switch sim.p.nameModel
    
    case 'chemostat'
        for iTime = 1:length(sim.t)
            u = [ sim.N(iTime) sim.DOC(iTime), sim.B(iTime,:) ];
            rates = getRates(p, u, sim.L, sim.T );
            % Losses from HTL:
            losses = losses + ...
                sum((1-reminHTL)*rates.mortHTL.*squeeze(sim.B(iTime,:)'))/rhoCN*dt(iTime);
            % Gains from diffusion:
            gains = gains + p.d*(p.u0(1)-sim.N(iTime)) * dt(iTime);
        end
        %
        % Calculate total budget:
        %
        accumulation = (sim.N(end)+sum(sim.B(end,:)/rhoCN)) - ...
            (sim.N(1)+sum(sim.B(1,:)/rhoCN));
        dNdt = (accumulation - gains + losses)/1000*p.depthProductiveLayer/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = (accumulation - gains + losses) / sim.N(end)/sim.t(end)*365; % Fraction per year
        
    case 'watercolumn'
        for iTime = 1:length(sim.t)
            for iDepth = 1:length(sim.z)
                u = [sim.N(iDepth,iTime) sim.DOC(iDepth,iTime), ...
                    sim.B(iDepth,:,iTime) ];
                rates = getRates(p, u, sim.L(iDepth,iTime), sim.T(iDepth,iTime) );
                % Losses from HTL:
                losses = losses + ...
                    sum((1-reminHTL)*rates.mortHTL.*squeeze(sim.B(iDepth,:,iTime))')/rhoCN*dt(iTime) ...
                    * sim.dznom(iDepth)/1000; %  gN/m2/day
            end
        end
        %
        % Calculate total budget:
        %
        accumulation = ...
            (sim.N(:,end)+sum(sim.B(:,:,end),2)/rhoCN) - ...
            (sim.N(:,1)+sum(sim.B(:,:,1),2)/rhoCN);
        accumulation = sum( accumulation.*sim.dznom )/1000;
        
        dNdt = (accumulation + losses)/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = (accumulation + losses) / (sum(sim.N(:,end).*sim.dznom)/1000);
        
    case 'global'
        losses = 0*sim.t;
        rates = getRates(p, p.u0, 100, 10 ); % Assume that HTL mortality does not vary!!
        load(sim.p.pathGrid,'dv');
        for iTime = 1:length(sim.t)
            for iDepth = 1:length(sim.z)
                for j = 1:length(sim.x)
                    for k = 1:length(sim.y)
                        if ~isnan(sim.N(j,k,iDepth,iTime))
                            %u = [sim.N(j,k,iDepth,iTime), sim.DOC(j,k,iDepth,iTime), ...
                            %    squeeze(sim.B(j,k,iDepth,:,iTime))' ];
                            %rates = getRates(p, u, sim.L(j,k,iDepth,iTime), sim.T(j,k,iDepth,iTime) );
                            % Losses from HTL:
                            losses(iTime) = losses(iTime) + ...
                                sum((1-reminHTL)*rates.mortHTL.*squeeze(sim.B(j,k,iDepth,:,iTime)))/rhoCN ...
                                * dv(j,k,iDepth)/1000 * dt(iTime); %  gN/day
                        end  
                    end
                end
            end
        end
        losses = sum(losses);
        %
        % Calculate total budget:
        %
        totN_0 = 0;
        totN_end = 0;
        for iDepth = 1:length(sim.z)
            for j = 1:length(sim.x)
                for k = 1:length(sim.y)
                    if ~isnan(sim.N(j,k,iDepth,iTime))
                        totN_0 = totN_0 + (sim.N(j,k,iDepth,1)  +sum(sim.B(j,k,iDepth,:,1),4)/rhoCN) * dv(j,k,iDepth)/1000;
                        totN_end = totN_end + (sim.N(j,k,iDepth,end)  +sum(sim.B(j,k,iDepth,:,end),4)/rhoCN) * dv(j,k,iDepth)/1000;                        
                    end
                end
            end
        end
        accumulation = sum( totN_end - totN_0 );
        
        dNdt = (accumulation + losses)/3.6e14/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = (accumulation + losses) / totN_end;
        
end
%
% Print result:
%
if bVerbose
    fprintf("N balance:\n");
    fprintf("  absolute: %f gN/m2/yr\n", dNdt);
    fprintf("  relative: %f 1/yr\n", dNdt_per_N);
end