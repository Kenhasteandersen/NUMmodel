%
% Plot Sheldon biomass spectrum. The biomasses are normalised by the width
% of the cell and multiplied by mass.
%
function panelSpectrum(sim, ixTime)

arguments
    sim struct;
    ixTime {mustBeInteger} = length(sim.t);
end
p = sim.p;

sLegend = {};

mMin = 100;
mMax = 0;

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    mMin = min([mMin, m]);
    mMax = max([mMax, m]);
    
    loglog(m, sim.B(ixTime,ix-p.idxB+1)./p.mDelta(ix).*m, 'linewidth',2)
    hold on
    
    sLegend{iGroup} = p.nameGroup{iGroup};
end
ylim([0.0001,500])
xlim([mMin, mMax])
hold off

xlabel('Mass ({\mu}gC)')
ylabel('Sheldon biomass ({\mu}gC/L)')

legend(sLegend, 'location','southeast','box','off')