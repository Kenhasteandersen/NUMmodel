
function panelTrophicLevel(p,rates)

lambda=calcTrophicLevel(p,rates);
nameGroup=p.nameGroup;

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);
    set(gca, 'XScale', 'log')
    hold on
    %
    % Generalists :
    %
    if ((p.typeGroups(iGroup)==1) || (p.typeGroups(iGroup)==5))
        loglog(m,lambda(ix),'linewidth',2,'Color',p.colGroup{iGroup})  
    %
    % Diatoms:
    %
    elseif ((p.typeGroups(iGroup)==3) || (p.typeGroups(iGroup)==4))
        loglog(m,lambda(ix),'linewidth',2,'Color',p.colGroup{iGroup})
    %
    % Copepods:
    %
    elseif ((p.typeGroups(iGroup)==10) || (p.typeGroups(iGroup)==11))
        loglog(m,lambda(ix),'linewidth',2,'Color',p.colGroup{iGroup})
    
    elseif p.typeGroups(iGroup)==100
        nameGroup(iGroup)=[];
    end
    
end

hold off

ylim([0.9 max(lambda)+0.3])
xlim(calcXlim(p))
ylabel('Trophic Level')
xlabel('Mass ({\mu}gC)')
legend(nameGroup,'Location','northwest','box','off')
