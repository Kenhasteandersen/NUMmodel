function panelTime(sim)
t = sim.t;

semilogy(t, sim.N,'b-')
hold on
semilogy(t, sim.DOC,'color',[181 100 30]/256)
textLegend = {'N','DOC'};
if isfield(sim,'Si')
    semilogy(t, sim.Si,'color',[181 180 0]/256)
    textLegend{3} = 'Si';
end

for iGroup = 1:sim.p.nGroups
    semilogy(t, sim.Bgroup(:,iGroup),'k-', 'linewidth',0.5*iGroup);
    textLegend{iGroup+sim.p.idxB-1} = sprintf('%2i',iGroup);
end
ylim([0.1, 10000])
axis('tight')

legend(textLegend,'location','eastoutside')

ylabel('Concentrations')
xlabel('Time (days)')