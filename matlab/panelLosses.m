function panelLosses(p,rates)

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    semilogx(m, rates.mortpred(ix), 'r-','linewidth',2)
    hold on
    semilogx(m, p.Jresp(ix)./m, 'k-', 'linewidth',2)
%    loglog(m, rates.mortStarve(ix), 'b-o','linewidth',2)
end
loglog(p.m, p.mortHTLm, 'm-','linewidth',2)

hold off
xlim([min(p.m), max(p.m)])

ylabel('Losses (day^{-1})')
xlabel('Mass ({\mu}gC)')
legend({'Predation','Respiration','HTL'}, ...
    'location','northwest','box','off')