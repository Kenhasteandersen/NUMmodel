function panelTime(sim)
t = sim.t;

semilogy(t, sim.N,'b-')
hold on
semilogy(t, sim.DOC,'color',[181 100 30]/256)
if isfield(sim,'Si')
    semilogy(t, sim.DOC,'color',[181 100 30]/256)
end

textLegend = {'N','DOC'};
for iGroup = 1:sim.p.nGroups
    semilogy(t, sim.Bgroup(:,iGroup),'k-', 'linewidth',0.5*iGroup);
    textLegend{iGroup+2} = sprintf('%2i',iGroup);
end
ylim([0.1, 10000])
axis('tight')

legend(textLegend,'location','eastoutside')

ylabel('Concentrations')
xlabel('Time (days)')