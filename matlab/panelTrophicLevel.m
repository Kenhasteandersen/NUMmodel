function [lambda, lambdaHTL] = panelTrophicLevel(p, B, rates, ax)

arguments
    p struct;
    B;
    rates struct;
    ax = [];
end

% Use gca if no axis is given
if isempty(ax)
    ax = gca;
end

% Compute trophic levels
[lambda, lambdaHTL] = calcTrophicLevel(p, B, rates);

% Prepare legend
name = {};
presence = zeros(1, 11);

% Clear and set up axis
cla(ax);
hold(ax, 'on');

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB + 1;
    m = p.m(ix + p.idxB - 1);
    
    loglog(ax, m, lambda(ix), 'linewidth', 2, 'Color', p.colGroup{iGroup});
    
    % Legend entries
    if (p.typeGroups(iGroup) == 1 || p.typeGroups(iGroup) == 5) && presence(p.typeGroups(iGroup)) == 0 
        name{end+1} = 'Generalists';
        presence(p.typeGroups(iGroup)) = 1;

    elseif (p.typeGroups(iGroup) == 3 || p.typeGroups(iGroup) == 4) && presence(p.typeGroups(iGroup)) == 0 
        name{end+1} = 'Diatoms';
        presence(p.typeGroups(iGroup)) = 1;

    elseif p.typeGroups(iGroup) == 10 && presence(p.typeGroups(iGroup)) == 0 
        name{end+1} = 'Passive copepod';
        presence(p.typeGroups(iGroup)) = 1;

    elseif p.typeGroups(iGroup) == 11 && presence(p.typeGroups(iGroup)) == 0 
        name{end+1} = 'Active copepod';
        presence(p.typeGroups(iGroup)) = 1;
    end
end

hold(ax, 'off');

ylim(ax, [0.9, max(lambda) + 0.3]);
xlim(ax, calcXlim(p));
ylabel(ax, 'Trophic Level');
xlabel(ax, 'Mass ({\mu}g_C)');
set(ax, 'XScale', 'log'); % 👈 Enforce logarithmic scale

legend(ax, name, 'Location', 'eastoutside', 'box', 'off');

end
