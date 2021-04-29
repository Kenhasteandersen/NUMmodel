%
% Plot Sheldon biomass spectrum. The biomasses are normalised by the width
% of the cell and multiplied by mass.
%
function panelSpectrum(sim)
p = sim.p;
ixTime = length(sim.t);

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    
    loglog(m, sim.B(ixTime,ix-p.idxB+1)./p.mDelta(ix).*m, 'linewidth',2)
    hold on
end
ylim([0.0001,500])
xlim([min(sim.p.m), max(sim.p.m)])
hold off

xlabel('Mass ({\mu}gC)')
ylabel('Sheldon biomass ({\mu}gC/L)')