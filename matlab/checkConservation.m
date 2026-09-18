%
% Calculates the nitrogen balance of the course of the simulation.
%
% Note: in the 'global' case all losses and the bottom-BC input are measured
% directly in simulateGlobal, so the residual is the non-conservation of the
% transport step. "Loss to HTL" is then the export from the biogeochemistry,
% i.e. HTL and POM losses when there is no POM group, and zero with POM.
%
% Note: for a watercolumn WITHOUT POM the balance is only accurate to about
% 0.1 %/yr. The HTL and POM losses leave the system directly in that case,
% and they are reconstructed from getRates once per save, whereas the model
% integrates them at the internal time step p.dt. That quadrature error is
% a few percent of the loss terms. With POM present the losses are measured
% directly and the balance closes to ~1e-6 /yr.
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
S = inputRead;
fracHTL_to_N = S.general.fracHTL_to_N;
rhoCN = S.general.rhoCN;

p = sim.p;
gains = 0;
lossHTL = 0;
dt = diff(sim.t);
dt = reshape(dt, length(dt),1);
dt = [dt; dt(end)];

switch sim.p.nameModel

    case 'chemostat'
        loss = 0;
        for iTime = 1:(length(sim.t)-1)
            if ~isfield(sim.p,'idxSi')
                u = [ sim.N(iTime) sim.DOC(iTime), sim.B(iTime,:) ];
            else
                u = [ sim.N(iTime) sim.DOC(iTime), sim.Si(end), sim.B(iTime,:) ];
            end
            rates = getRates(p, u, mean(sim.L), sim.T );
            ixUni = findIxUnicellular(sim.p);

            B = squeeze(0.5*sum(sim.B(iTime:iTime+1,:))); % Interpolate B
            if ~sum(ismember(p.typeGroups,100))
                %
                % If pom is not present:
                %
                
                % Losses from HTL:
                lossHTL = lossHTL + ...
                    sum((1-fracHTL_to_N)*rates.mortHTL.*B')/rhoCN*dt(iTime);
                % Losses to POM:
                loss = loss + sum(rates.jPOM.*B')/rhoCN*dt(iTime);%sum(rates.mort2*(1-remin2).*B')/rhoCN*dt(iTime);
            else
                %
                % If POM is present:
                %
                ixPOM = p.ixStart(p.ixPOM):p.ixEnd(p.ixPOM);
                loss = loss + p.velocity(ixPOM).*u(ixPOM)/p.widthProductiveLayer /rhoCN*dt(iTime);
                %ixPOM = p.ixStart(ixGroupPOM):p.ixEnd(ixGroupPOM);
                %           loss = loss + sim(p.velocity(ixPOM).*u(ixPOM))/rhoCN*dt(iTime);
            end

            % Losses from diffusion:
            loss = loss + sim.bUnicellularloss*p.d*sum(B(ixUni))/rhoCN * dt(iTime);

            % Gains from diffusion:
            gains = gains + p.d*(p.uDeep(p.idxN)-0.5*(sim.N(iTime)+sim.N(iTime+1))) * dt(iTime);
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
        % The loss and production terms are rates (gN/m2/day), so they are
        % weighted by the time between saves:
        Ntot = reshape(sim.Ntot, [], 1);
        accumulation = Ntot-Ntot(1) + p.tSave*( ...
            -cumsum(sim.Nprod) + cumsum(sim.NlossHTL) + cumsum(sim.Nloss) );

        dNdt = (accumulation(end)-accumulation(1))/sim.t(end)*365; %gN/m2/yr
        dNdt_per_N = dNdt/sim.Ntot(end);
        lossHTL = p.tSave*cumsum(sim.NlossHTL)/sim.t(end)*365; %gN/m2/yr
        lossHTL = lossHTL(end);
        lossHTL_per_N = lossHTL/sim.Ntot(end);

    case 'global'
        %
        % All terms are totals over the ocean (gN) measured in simulateGlobal
        % over each save interval. The fluxes recorded at the first save cover
        % the interval from t=0 and are dropped, so that the budget runs from
        % the first to the last saved state. Any residual is what the
        % transport step fails to conserve.
        %
        if length(sim.t) < 2
            error('checkConservation needs at least two saved time points (decrease p.tSave).');
        end
        areaOcean = 3.6e14; % m2
        T = sim.t(end) - sim.t(1); % days
        accumulation = (sim.Ntot(end)-sim.Ntot(1)) ...
            - sum(sim.Nprod(2:end)) + sum(sim.Nloss(2:end));

        dNdt = accumulation/areaOcean/T*365; %gN/m2/yr
        dNdt_per_N = accumulation/sim.Ntot(end)/T*365; % 1/yr
        lossHTL = sum(sim.NlossHTL(2:end))/areaOcean/T*365; %gN/m2/yr
        lossHTL_per_N = sum(sim.NlossHTL(2:end))/sim.Ntot(end)/T*365; % 1/yr

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