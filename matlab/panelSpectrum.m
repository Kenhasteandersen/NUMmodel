%
% Plot Sheldon biomass spectrum. The biomasses are normalised by the log of
% the ratio between upper and lower masses in each bin
%
function panelSpectrum(sim, ixTime)

arguments
    sim struct;
    ixTime {mustBeInteger} = length(sim.t);
end
p = sim.p;

sLegend = {};

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    
    loglog(m, sim.B(ixTime,ix-p.idxB+1)./log(p.mUpper(ix)./p.mLower(ix)), 'linewidth',2)
    hold on
    
    sLegend{iGroup} = p.nameGroup{iGroup};
end
ylim([0.1,500])
xlim(calcXlim(sim.p))
hold off

xlabel('Mass ({\mu}gC)')
ylabel('Sheldon biomass ({\mu}gC/L)')

legend(sLegend, 'location','northeast','box','off')