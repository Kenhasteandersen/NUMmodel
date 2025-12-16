%
% Makes a panel with the trophic levels of each spectrum group
%
% In:
%  p - a structure with parameters from a setup
%  rates - a structure with jDOC, jLreal and jF rates (from getRates)
%
% Out:
%  lambda - the trophic level vector of all size groups that are not
%          nutrients.
%  lambdaHTL - The trophic level of the HTL calculated as 1 plus the biomass
%           weighted mean of the HTL consumption.
%

function [lambda, lambdaHTL] = panelTrophicLevel(p,B,rates)

arguments
    p struct;
    B;
    rates struct;
end

[lambda, lambdaHTL] = calcTrophicLevel(p,B,rates);
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
xlabel('Mass ({\mu}g_C)')
legend(name,'Location','eastoutside','box','off')
