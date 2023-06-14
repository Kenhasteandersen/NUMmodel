
function panelTrophicLevel(sim,rates)

lambda=calcTrophicLevel(sim,rates);
p=sim.p;
name={};
presence=zeros(1,11);

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
    m = p.m(ix+p.idxB-1);
    set(gca, 'XScale', 'log')
    hold on
    loglog(m,lambda(ix),'linewidth',2,'Color',p.colGroup{iGroup})
        
    %
    % Legend :
    %

    % Generalists :
    if (p.typeGroups(iGroup)==1 || p.typeGroups(iGroup)==5) && presence(p.typeGroups(iGroup))==0 
           name{end+1}='Generalists';
           presence(p.typeGroups(iGroup))=1;
       
    %Diatoms:
    elseif (p.typeGroups(iGroup)==3 || p.typeGroups(iGroup)==4) && presence(p.typeGroups(iGroup))==0 
           name{end+1}='Diatoms';
           presence(p.typeGroups(iGroup))=1;
    % Copepods:
    elseif (p.typeGroups(iGroup)==10)&& presence(p.typeGroups(iGroup))==0 
           name{end+1}='Passive copepod';
           presence(p.typeGroups(iGroup))=1;

    elseif (p.typeGroups(iGroup)==11)&& presence(p.typeGroups(iGroup))==0 
           name{end+1}='Active copepod';
           presence(p.typeGroups(iGroup))=1;
    end
    
end

hold off

ylim([0.9 max(lambda)+0.3])
xlim(calcXlim(p))
ylabel('Trophic Level')
xlabel('Mass ({\mu}gC)')
legend(name,'Location','northwest','box','off')
