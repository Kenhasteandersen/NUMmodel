function panelGains(p,rates, bDrawStrategies)

arguments
    p, rates struct;
    bDrawStrategies = false;
end

% brackground depending on trophic strategies
if bDrawStrategies
    calcTrophicStrategy(p, rates);
end


for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);

    %
    % Generalists and copepods:
    %
    if ((p.typeGroups(iGroup)~=3) && (p.typeGroups(iGroup)~=4))
        semilogx(m, rates.jFreal(ix), 'r-', 'linewidth',2)
        hold on
        semilogx(m, rates.jN(ix), 'b-','linewidth',2)
        semilogx(m, rates.jLreal(ix), 'g-','linewidth',2)
        semilogx(m, rates.jDOC(ix), 'color',[181 100 30]/256,'linewidth',2)
        semilogx(m, rates.jMax(ix), 'k:')
        % Total:
        semilogx(m, rates.jTot(ix), 'k-', 'linewidth',2)

    end
    %
    % Diatoms:
    %
    if ((p.typeGroups(iGroup)==3) || (p.typeGroups(iGroup)==4))
        semilogx(m, rates.jSi(ix), 'color',[181 180 0]/256,'linewidth',2)
        semilogx(m, rates.jN(ix), 'b--','linewidth',2)
        semilogx(m, rates.jLreal(ix), 'g--','linewidth',2)
        semilogx(m, rates.jDOC(ix), 'color',[181 100 30]/256,'linewidth',2)
        semilogx(m, rates.jMax(ix), 'k:')
        semilogx(m, rates.jTot(ix), 'k--', 'linewidth',2)
    end
    %
    % Rates for copepods:
    %
    %if (p.typeGroups(iGroup) >= 10)
    %    semilogx(m, rates.jF(ix), 'r-', 'linewidth',2)
    %    hold on
    %    semilogx(m, p.jFmax(ix),'k:')
    %end

end



hold off
ylim([0 2])
xlim(calcXlim(p))
ylabel('Gains (day^{-1})')

if isfield(p, 'ixSi')
    cap={'Feeding','N','Light','Si','DOC','Max. growth rate','Growth rate'};
    legend(cap, ...
        'location','northwest','box','off')
else
    cap={'Feeding','N','Light','DOC','Max. growth rate','Growth rate'};
    legend(cap, ...
        'location','northwest','box','off')
end
