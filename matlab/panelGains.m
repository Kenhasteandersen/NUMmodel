function panelGains(p,rates)

for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    semilogx(m, rates.JF(ix)./m, 'r-o', 'linewidth',2)
    hold on
    % Rates for generalists:
    if (p.typeGroups(iGroup)==1)
       semilogx(m, rates.JN(ix)./m, 'b-','linewidth',2)
       semilogx(m, rates.JLreal(ix)./m, 'g-','linewidth',2)
       semilogx(m, rates.JDOC(ix)./m, 'color',[181 100 30]/256,'linewidth',2)
       semilogx(m, p.pGeneralists.Jmax./m, 'k:')
    end
    % Total:
    semilogx(m, rates.Jtot(ix)./m, 'k-', 'linewidth',2)
    % Rates for copepods:
    if (p.typeGroups(iGroup)==2)
        semilogx(m, p.JFmax(ix)./m,'k:')
    end
end

hold off
ylim([0 2])
xlim([min(p.m), max(p.m)])
ylabel('Gains (day^{-1})')

legend({'Feeding','N','Light','DOC','Max. growth rate','Growth rate'}, ...
    'location','northwest','box','off')
