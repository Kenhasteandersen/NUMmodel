function plot_diagnostics(p,sim,rates)


figure
subplot(3,1,1)
panelSpectrum(sim)


subplot(3,1,2)
for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    semilogx(m, rates.jF(ix)', 'r-', 'linewidth',2)
    hold on
    % Rates for generalists:
    if (p.typeGroups(iGroup)==1)
       semilogx(m, rates.jN(ix)', 'b-','linewidth',2)
       semilogx(m, rates.jL(ix)', 'g-','linewidth',2)
%        semilogx(m, rates.jDOC(ix)./m, 'color',[181 100 30]/256,'linewidth',2)
%        semilogx(m, p.pGeneralists.jTot./m, 'k:')
    end
    % Total:
    semilogx(m, rates.jTot(ix)', 'k-', 'linewidth',2)
    % Rates for copepods:
    if (p.typeGroups(iGroup)==2)
%         semilogx(m, p.JFmax(ix)./m,'k:')
    end
end

hold off
ylim([0 2])
xlim([min(p.m), max(p.m)])
ylabel('Gains (day^{-1})')

legend({'Feeding','N','Light','Growth rate'}, ...
    'location','northwest','box','off')



subplot(3,1,3)
for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    semilogx(m, rates.mortpred(ix)', 'r-','linewidth',2)
    hold on
    semilogx(m, rates.mortHTL(ix)', 'k-', 'linewidth',2)
%     semilogx(m, rates.mort2, 'b-','linewidth',2)
%    loglog(m, rates.mortStarve(ix), 'b-o','linewidth',2)
end
% loglog(p.m, p.mortHTLm, 'm-','linewidth',2)

hold off
xlim([min(p.m), max(p.m)])

ylabel('Losses (day^{-1})')
xlabel('Mass ({\mu}gC)')
legend({'Predation','HTL'}, ...
    'location','northwest','box','off')
