function panelLosses(p, rates, bDrawStrategies, ax)
arguments
    p, rates struct;
    bDrawStrategies = false;
    ax = [];
end

% Use gca if no axis is provided
if isempty(ax)
    ax = gca;
end

cla(ax);
hold(ax, 'on');

% Background depending on trophic strategies
if bDrawStrategies
    calcTrophicStrategy(p, rates,ax);
end

for iGroup = 1:p.nGroups
    if ((p.typeGroups(iGroup) == 3) || (p.typeGroups(iGroup) == 4))
        sLinetype = "--"; % For diatoms
    else
        sLinetype = "-"; % For all others
    end

    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB + 1;
    m = p.m(ix + p.idxB - 1);

    semilogx(ax, m, rates.mortpred(ix), 'r', 'linewidth', 2, 'linestyle', sLinetype);
    semilogx(ax, m, rates.jRespTot(ix), 'k', 'linewidth', 2, 'linestyle', sLinetype);
    semilogx(ax, m, rates.mort2(ix), 'b', 'linewidth', 2, 'linestyle', sLinetype);
    semilogx(ax, m, rates.jLossPassive(ix), 'color', [0 0.5 0], 'linewidth', 2, 'linestyle', sLinetype);
    loglog(ax, m, rates.mortHTL(ix), 'm', 'linewidth', 2, 'linestyle', sLinetype);
end

% Show losses from chemostat
if strcmp(p.nameModel, 'chemostat') && (isnan(p.seasonalOptions.lat_lon) && p.seasonalOptions.seasonalAmplitude == 0)
    semilogx(ax, m, p.d * m ./ m, '--');
end

hold(ax, 'off');
xlim(ax, calcXlim(p));
ylabel(ax, 'Losses (day^{-1})');
xlabel(ax, 'Mass ({\mu}g_C)');
set(ax, 'XScale', 'log'); % Ensures logarithmic scale

legend(ax, {'Predation', 'Respiration', 'Viral lysis', 'Passive exudation', 'HTL'}, ...
    'location', 'eastoutside', 'box', 'off');
end
