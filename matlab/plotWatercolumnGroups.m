%
% Plot a vertical profile of the total biomass of each group in a
% simulation.  Currently very basic. Only works for a watercolumn
% simulation.
%
% In:
%  sim - the simulation structure
%  tDay - The day for which to plot the profile
%
function plotWatercolumnGroups(sim, tDay)

p = sim.p;

if ~strcmp(p.nameModel,'watercolumn')
    stop('Only works for water column models.')
else
    clf
    iDay = find(abs(tDay-sim.t) == min(abs(tDay-sim.t)));

    z = sim.z;

    plot(sim.N(iDay,:), -z, 'b-','linewidth',2)
    hold on

    if isfield(sim,'Si')
        plot(sim.Si(iDay,:),-z,'y-','linewidth',2)
    end

    for iGroup = 1:p.nGroups
        ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
        plot(squeeze(sum(sim.B(iDay,:,ix),3)), -z, "Color",p.colGroup{iGroup},'linewidth',2)
    end

    ylabel('Depth (m)')
    xlabel('Concentration ({\mu}g/l)')

    xlim([0 200])
    ylim([-max(z) 0])
end