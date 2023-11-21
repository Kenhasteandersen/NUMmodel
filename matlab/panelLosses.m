function panelLosses(p,rates, bDrawStrategies)
arguments
    p, rates struct;
    bDrawStrategies = false;
end

% background depending on trophic strategies
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
    if ((p.typeGroups(iGroup)==3) || (p.typeGroups(iGroup)==4))
        sLinetype = "--"; % For diatoms
    else
        sLinetype = "-"; % For all others
    end

    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);

    nexttile(1)
    semilogx(m, rates.mortpred(ix), 'color',p.colGroup{iGroup},'linewidth',2, 'linestyle',sLinetype)
    hold on
    sLegend{iGroup,1} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
    %semilogx(m, rates.jR(ix), 'k-.', 'linewidth',2)
    nexttile(2)
    semilogx(m, rates.jRespTot(ix), 'color',p.colGroup{iGroup}, 'linewidth',2, 'linestyle',sLinetype)
    hold on
    sLegend{iGroup,2} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
    nexttile(3)
    semilogx(m, rates.mort2(ix), 'color',p.colGroup{iGroup},'linewidth',2, 'linestyle',sLinetype)
    hold on
    sLegend{iGroup,3} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
    nexttile(4)
    semilogx(m, rates.jLossPassive(ix), 'color',p.colGroup{iGroup},'linewidth',2, 'linestyle',sLinetype)
    hold on
    sLegend{iGroup,4} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
    %    loglog(m, rates.mortStarve(ix), 'b-o','linewidth',2)
    nexttile(5)
    loglog(m, rates.mortHTL(ix), 'color',p.colGroup{iGroup},'linewidth',2, 'linestyle',sLinetype)
    hold on
    sLegend{iGroup,5} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
end
%
% Show losses from chemostat
%
if strcmp(p.nameModel,'chemostat') & (isnan(p.seasonalOptions.lat_lon) & p.seasonalOptions.seasonalAmplitude==0)
    nexttile(6)
    semilogx(m, p.d*m./m,'--','color',p.colGroup{iGroup})
    hold on
    sLegend{iGroup,6} =[p.nameGroup{iGroup},' pal=',num2str(palatability(iGroup)),' alphaF=',num2str(alphaF(iGroup))];
end

tit={'Predation','Respiration','Viral lysis','Passive exudation','HTL','chemostat'};
for itnr=1:6
    nexttile(itnr);hold on
    % ylim([0 2])
    xlim(calcXlim(p))
    ylabel('Losses (day^{-1})')
     xlabel('Mass ({\mu}gC)')
    title(tit{itnr})
    if itnr==5
        legend(sLegend(:,5), ...
            'location','best','box','off')
    end
end

%hold off
%xlim(calcXlim(p))

% ylabel('Losses (day^{-1})')
% xlabel('Mass ({\mu}gC)')
% legend({'Predation','Respiration','Viral lysis','Passive exudation','HTL'}, ...
%     'location','northwest','box','off')