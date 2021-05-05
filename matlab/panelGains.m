function panelGains(p,rates)

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);
    semilogx(m, rates.jF(ix), 'r-', 'linewidth',2)
    hold on
    % Rates for unicellulars:
    if (p.typeGroups(iGroup)<10)
       semilogx(m, rates.jN(ix), 'b-','linewidth',2)
       semilogx(m, rates.jLreal(ix), 'g-','linewidth',2)
       semilogx(m, rates.jDOC(ix), 'color',[181 100 30]/256,'linewidth',2)
       semilogx(m, rates.jMax(ix), 'k:')
    end
    
    % Total:
    semilogx(m, rates.jTot(ix), 'k-', 'linewidth',2)
    % Rates for copepods:
    %if (p.typeGroups(iGroup)~=2)
    %    semilogx(m, p.jFmax(ix),'k:')
    %end
end

hold off
ylim([0 2])
xlim([min(p.m), max(p.m)])
ylabel('Gains (day^{-1})')

legend({'Feeding','N','Light','DOC','Max. growth rate','Growth rate'}, ...
    'location','northwest','box','off')
