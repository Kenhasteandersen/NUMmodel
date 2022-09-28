function panelLosses(p,rates)

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);
    semilogx(m, rates.mortpred(ix), 'r-','linewidth',2)
    hold on
    semilogx(m, rates.jR(ix), 'k-', 'linewidth',2)
    semilogx(m, rates.mort2(ix), 'b-','linewidth',2)
    semilogx(m, rates.jLossPassive(ix), 'color',[0 0.5 0],'linewidth',2)
%    loglog(m, rates.mortStarve(ix), 'b-o','linewidth',2)
    loglog(m, rates.mortHTL(ix), 'm-','linewidth',2)
end
%
% Show losses from chemostat
%
if strcmp(p.nameModel,'chemostat') & (isnan(p.seasonalOptions.lat_lon) & p.seasonalOptions.seasonalAmplitude==0)
    semilogx(m, p.d*m./m,'--')
end

hold off
xlim(calcXlim(p))

ylabel('Losses (day^{-1})')
xlabel('Mass ({\mu}gC)')
legend({'Predation','Respiration','Viral lysis','Passive exudation','HTL'}, ...
    'location','northwest','box','off')