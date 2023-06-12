function panelLosses(p,rates)

% brackground depending on trophic strategies

[strategy, col] = calcTrophicStrategy(rates);

Ngroup=1; %select the group for which the background is drawn
ix = (p.ixStart(Ngroup):p.ixEnd(Ngroup));
m = p.m(ix);
color=col(ix-(p.idxB-1),:);
colori=color(1,:);
Xlim=calcXlim(p);
Xmin=Xlim(1);
Xmax=Xlim(2);
rectangle(Position=[Xmin,0,(m(2)+m(1))/2-Xmin,2], FaceColor=colori, EdgeColor=colori);
hold on
for i = 2:length(ix)-1
    colori=color(i,:);
    step=(m(i)+m(i+1))/2-(m(i-1)+m(i))/2;
    rectangle(Position=[(m(i-1)+m(i))/2,0,step,2], FaceColor=colori, EdgeColor=colori);
end
colori=color(length(ix),:);
rectangle(Position=[(m(length(ix)-1)+m(length(ix)))/2,0,Xmax-(m(length(ix)-1)+m(length(ix)))/2,2], FaceColor=colori, EdgeColor=colori);


set(gca,'XScale', 'log', "Layer", "top")


for iGroup = 1:p.nGroups
    if ((p.typeGroups(iGroup)==3) || (p.typeGroups(iGroup)==4))
        sLinetype = "--"; % For diatoms
    else
        sLinetype = "-"; % For all others
    end

    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);

    semilogx(m, rates.mortpred(ix), 'r','linewidth',2, 'linestyle',sLinetype)
    hold on
    %semilogx(m, rates.jR(ix), 'k-.', 'linewidth',2)
    semilogx(m, rates.jRespTot(ix), 'k', 'linewidth',2, 'linestyle',sLinetype)
    semilogx(m, rates.mort2(ix), 'b','linewidth',2, 'linestyle',sLinetype)
    semilogx(m, rates.jLossPassive(ix), 'color',[0 0.5 0],'linewidth',2, 'linestyle',sLinetype)
    %    loglog(m, rates.mortStarve(ix), 'b-o','linewidth',2)
    loglog(m, rates.mortHTL(ix), 'm','linewidth',2, 'linestyle',sLinetype)
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