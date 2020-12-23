function panelSpectrum(sim)
p = sim.p;
ixTime = length(sim.t);

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    
    loglog(p.m(ix), sim.B(ixTime,ix-2), 'linewidth',2)
    hold on
end
ylim([0.001,500])
hold off

xlabel('Mass ({\mu}gC)')
ylabel('Biomass ({\mu}gC/L)')