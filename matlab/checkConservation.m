%
% Calculates the nitrogen balance of the course of the simulation.
%
% Note: in the 'chemostat' case the budget terms are integrated by the ODE
% solver together with the state, so the balance closes to the solver
% tolerance (~1e-8 /yr with the default tolerances).
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
rhoCN = S.general.rhoCN;

p = sim.p;

switch sim.p.nameModel

    case 'chemostat'
        %
        % Nprod, Nloss and NlossHTL are integrated along with the state in
        % simulateChemostat (cumulative from t=0, mugN/l), so the budget
        % closes to the tolerance of the ODE solver.
        %
        Ntot = sim.N + sum(sim.B,2)/rhoCN;
        accumulation = (Ntot(end)-Ntot(1)) - sim.Nprod(end) + sim.Nloss(end);
        T = sim.t(end);

        dNdt = accumulation/1000*p.widthProductiveLayer/T*365; %gN/m2/yr
        dNdt_per_N = accumulation/Ntot(end)/T*365; % 1/yr
        lossHTL = sim.NlossHTL(end)/1000*p.widthProductiveLayer/T*365; %gN/m2/yr
        lossHTL_per_N = sim.NlossHTL(end)/Ntot(end)/T*365; % 1/yr

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