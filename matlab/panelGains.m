function h = panelGains(p, rates, bDrawStrategies, ax)
arguments
    p struct
    rates struct
    bDrawStrategies = false
    ax = []
end

% Activate axis
if isempty(ax)
    ax = gca;
end
axes(ax); % Activate

cla(ax); % Clear axis content
hold(ax, 'on');

% Output handle storage
h.lines = [];
h.legend = [];

% Background strategies
if bDrawStrategies
    calcTrophicStrategy(p, rates); % No handle returned here
end

% Plot gains per group
for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB + 1;
    m = p.m(ix + p.idxB - 1);

    if ((p.typeGroups(iGroup) ~= 3) && (p.typeGroups(iGroup) ~= 4))
        h.lines(end+1) = semilogx(ax, m, rates.jFreal(ix), 'r-', 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jN(ix), 'b-', 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jLreal(ix), 'g-', 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jDOC(ix), 'Color', [181 100 30]/256, 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jMax(ix), 'k:');
        h.lines(end+1) = semilogx(ax, m, rates.jTot(ix), 'k-', 'LineWidth', 2);
    end

    if ((p.typeGroups(iGroup) == 3) || (p.typeGroups(iGroup) == 4))
        h.lines(end+1) = semilogx(ax, m, rates.jSi(ix), 'Color', [181 180 0]/256, 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jN(ix), 'b--', 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jLreal(ix), 'g--', 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jDOC(ix), 'Color', [181 100 30]/256, 'LineWidth', 2);
        h.lines(end+1) = semilogx(ax, m, rates.jMax(ix), 'k:');
        h.lines(end+1) = semilogx(ax, m, rates.jTot(ix), 'k--', 'LineWidth', 2);
    end
end

hold(ax, 'off');

ylim(ax, [0 2]);
xlim(ax, calcXlim(p));
ylabel(ax, 'Gains (day^{-1})');

if isfield(p, 'ixSi')
    cap = {'Feeding','N','Light','Si','DOC','Max. growth rate','Growth rate'};
else
    cap = {'Feeding','N','Light','DOC','Max. growth rate','Growth rate'};
end

h.legend = legend(ax, cap, 'Location', 'eastoutside', 'Box', 'off');
end
