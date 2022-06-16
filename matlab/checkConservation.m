%
% Calculates the nitrogen balance of the course of the simulation.
% Note: ONLY WORKS WITH CONSTANT HTL MORTALITY (NOT "QUADRATIC")
%
% In:
%  sim: simulation structure
%  bVerbose: whether to print out the balance on the terminal (default true).
%
function [dNdt, dNdt_per_N] = checkConservation(sim, bVerbose)

arguments
    sim struct;
    bVerbose logical = true; % Whether to print out the solution
end
%
% Constants:
%
fracHTL_to_N = 0.5;
rhoCN = 5.68;
remin2 = 0.5;

p = sim.p;
gains = 0;
lossHTL = 0;
dt = diff(sim.t);

switch sim.p.nameModel

    case 'chemostat'
        loss = 0;
        for iTime = 1:(length(sim.t)-1)
            u = [ sim.N(iTime) sim.DOC(iTime), sim.B(iTime,:) ];
            rates = getRates(p, u, mean(sim.L), sim.T );
            ixUni = findIxUnicellular(sim.p);
            if ~sum(ismember(p.typeGroups,100))
                B = squeeze(0.5*sum(sim.B(iTime:iTime+1,:))); % Interpolate B
                % Losses from HTL:
                lossHTL = lossHTL + ...
                    sum((1-fracHTL_to_N)*rates.mortHTL.*B')/rhoCN*dt(iTime);
                % Losses to POM:
                loss = loss + sum(rates.jPOM.*B')/rhoCN*dt(iTime);%sum(rates.mort2*(1-remin2).*B')/rhoCN*dt(iTime);
            else
                ixPOM = p.ixStart(ixGroupPOM):p.ixEnd(ixGroupPOM);
                loss = loss + sim(p.velocity(ixPOM).*u(ixPOM))/rhoCN*dt(iTime);
            end

            % Losses from diffusion:
            loss = loss + sim.bUnicellularloss*p.d*sum(B(ixUni))/rhoCN * dt(iTime);

            % Gains from diffusion:
            gains = gains + p.d*(p.uDeep-0.5*(sim.N(iTime)+sim.N(iTime+1))) * dt(iTime);
        end
        %
        % Calculate total budget:
        %
        accumulation = (sim.N(end)+sum(sim.B(end,:)/rhoCN)) - ...
            (sim.N(1)+sum(sim.B(1,:)/rhoCN));
        dNdt = (accumulation - gains + lossHTL + loss)/1000*p.widthProductiveLayer/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = (accumulation - gains + lossHTL+loss) / sim.N(end)/sim.t(end)*365; % Fraction per year
        lossHTL = lossHTL/1000*p.widthProductiveLayer/sim.t(end)*365;
        lossHTL_per_N = lossHTL/sim.N(end);

    case 'watercolumn'
        %
        % Calculate total budget:
        %
        accumulation = sim.Ntot-sim.Ntot(1) - cumsum(sim.Nprod) + cumsum(sim.NlossHTL) + cumsum(sim.Nloss);

        dNdt = (accumulation(end)-accumulation(1))/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = dNdt/sim.Ntot(end);
        lossHTL = cumsum(sim.NlossHTL)/sim.t(end)*365; %gN/m2/yr
        lossHTL = lossHTL(end);
        lossHTL_per_N = lossHTL/sim.Ntot(end);

    case 'global'
        losses = 0*sim.t;
        lossHTL = losses;
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
                            lossHTL(iTime) = lossHTL(iTime) + ...
                                sum((1-fracHTL_to_N)*rates.mortHTL.*squeeze(sim.B(j,k,iDepth,:,iTime)))/rhoCN ...
                                * dv(j,k,iDepth)/1000 * dt(iTime); %  gN/day
                            % Quadratic losses:
                            losses(iTime) = losses(iTime) + ...
                                (1-remin2)*sum(rates.mort2.*squeeze(sim.B(j,k,iDepth,:,iTime)))/rhoCN ...
                                * dv(j,k,iDepth)/1000 * dt(iTime); %  gN/day
                        end
                    end
                end
            end
        end
        losses = sum(losses);
        lossHTL = sum(lossHTL);
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

        dNdt = (accumulation + losses + lossHTL)/3.6e14/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = (accumulation + losses + lossHTL) / totN_end;
        lossHTL = lossHTL/3.6e14/sim.t(end)*365; %gN/m2/yr
        lossHTL_per_N = lossHTL/sim.Ntot(end);

end
%
% Print result:
%
if bVerbose
    fprintf("N balance:\n");
    fprintf("  loss to HTL: %f gN/m2/yr\n", lossHTL);
    fprintf("  relative HTL loss: %f 1/yr\n", lossHTL_per_N);
    fprintf("  absolute balance: %f gN/m2/yr\n", dNdt);
    fprintf("  relative balance: %f 1/yr\n", dNdt_per_N);
end