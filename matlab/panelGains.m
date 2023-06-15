function panelGains(p,rates)

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
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);
    
    %
    % Generalists and copepods:
    %
    if ((p.typeGroups(iGroup)~=3) && (p.typeGroups(iGroup)~=4))
       semilogx(m, rates.jF(ix), 'r-', 'linewidth',2)
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
