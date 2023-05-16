%
% Plots panels with total biomass of each group (g/m2) from a global simulation
%
% In:
%  sim: the simulation to plot
%  iTime: (not needed) the time step to plot (defaults to the last).
%  sProjection: default 'fast'.
%  bAverageTime=false: whether to average over time (defaults to false).
%
function plotGlobalGroups(sim, iTime, options)
arguments
    sim struct;
    iTime = length(sim.t);
    options.sProjection string = 'fast';
    options.bAverageTime = false;
end

clf
tiledlayout(sim.p.nGroups,1)

for i = 1:sim.p.nGroups
    nexttile
    % Extract the group:
    ix = (sim.p.ixStart(i):sim.p.ixEnd(i)) -sim.p.idxB+1;
    B = calcIntegrateGlobal(sim, sim.B(:,:,:,:,ix), options.bAverageTime);
    % Average over time:
    if ~options.bAverageTime
        B = B(iTime,:,:);
    end
    % Plot the group:
    c = panelGlobal(sim.x,sim.y, log10(B),[-1 2], ...
        sTitle=sim.p.nameGroup{i},...
        sProjection=options.sProjection);
    c.Label.String  = 'log10(g/m^2)';
end