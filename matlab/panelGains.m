function panelGains(p,rates, bDrawStrategies)

arguments
    p, rates struct;
    bDrawStrategies = false;
end

% brackground depending on trophic strategies
if bDrawStrategies
    calcTrophicStrategy(p, rates);
end
palatability1 = search_namelist('../input/input_generalists_simpleX.h','prokaryote','palatability');
alphaF1 = search_namelist('../input/input_generalists_simpleX.h','prokaryote','alphaF');
palatability = search_namelist('../input/input_generalists_simpleX.h','generalists_simple','palatability');
alphaF = search_namelist('../input/input_generalists_simpleX.h','generalists_simple','alphaF');
palatability=[palatability1,palatability];
alphaF=[alphaF1,alphaF];

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);

    %
    % Generalists and copepods:
    %
    if ((p.typeGroups(iGroup)~=3) && (p.typeGroups(iGroup)~=4))
        nexttile(4)
        semilogx(m, rates.jF(ix), '-', 'linewidth',2, 'color',p.colGroup{iGroup})
        hold on
        sLegend{iGroup,4} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
        nexttile(5)
        semilogx(m, rates.jN(ix), '-','linewidth',2, 'color',p.colGroup{iGroup})
        hold on
        sLegend{iGroup,5} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
        nexttile(3)
        semilogx(m, rates.jLreal(ix), '-','linewidth',2, 'color',p.colGroup{iGroup});
        hold on
        sLegend{iGroup,3} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
        nexttile(2)
        semilogx(m, rates.jDOC(ix),'color',p.colGroup{iGroup},'linewidth',2)
        hold on
        sLegend{iGroup,2} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
        nexttile(1)
        semilogx(m, rates.jMax(ix), ':','color',p.colGroup{iGroup})
        hold on
        % Total:
        semilogx(m, rates.jTot(ix), '-','color',p.colGroup{iGroup},'linewidth',2)
        sLegend{iGroup,1} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];

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



% hold off
% ylim([0 2])
% xlim(calcXlim(p))
% ylabel('Gains (day^{-1})')

tit={'Growth rate','DOC','Light','Feeding','N'};
for itnr=1:5
    nexttile(itnr);hold on
    % ylim([0 2])
    xlim(calcXlim(p))
    ylabel('Gains (day^{-1})')
    title(tit{itnr})
    if itnr==5
        legend(sLegend(:,5), ...
            'location','best','box','off')
    end
end



% if isfield(p, 'ixSi')
%     cap={'Feeding','N','Light','Si','DOC','Max. growth rate','Growth rate'};
%     legend(cap, ...
%         'location','northwest','box','off')
% else
%     cap={'Feeding','N','Light','DOC','Max. growth rate','Growth rate'};
%     legend(cap, ...
%         'location','northwest','box','off')
% end
