function lh = panelGainsTest(p, rates, bDrawStrategies, ax)
arguments
    p struct
    rates struct
    bDrawStrategies = false
    ax = []
end

% If no axis is provided, use the current one
if isempty(ax)
    ax = gca;
end

cla(ax); % Clear the axis
hold(ax, 'on');

% Initialize output handles
lh.lines = [];
lh.legend = [];

% Background depending on trophic strategies
if bDrawStrategies
    calcTrophicStrategy(p, rates); % Assumes it manages its own axis
end

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB + 1;
    m = p.m(ix + p.idxB - 1);

    % Generalists and copepods:
    if ((p.typeGroups(iGroup) ~= 3) && (p.typeGroups(iGroup) ~= 4))
        lh.lines(end+1) = semilogx(ax, m, rates.jFreal(ix), 'r-', 'linewidth', 2);            % Feeding
        lh.lines(end+1) = semilogx(ax, m, rates.jN(ix), 'b-', 'linewidth', 2);               % Nitrogen
        lh.lines(end+1) = semilogx(ax, m, rates.jLreal(ix), 'g-', 'linewidth', 2);           % Light
        lh.lines(end+1) = semilogx(ax, m, rates.jDOC(ix), 'color', [181 100 30]/256, 'linewidth', 2); % DOC
        lh.lines(end+1) = semilogx(ax, m, rates.jMax(ix), 'k:');                             % Max growth rate
        lh.lines(end+1) = semilogx(ax, m, rates.jTot(ix), 'k-', 'linewidth', 2);             % Total growth
    end

    % Diatoms:
    if ((p.typeGroups(iGroup) == 3) || (p.typeGroups(iGroup) == 4))
        lh.lines(end+1) = semilogx(ax, m, rates.jSi(ix), 'color', [181 180 0]/256, 'linewidth', 2);   % Silicate
        lh.lines(end+1) = semilogx(ax, m, rates.jN(ix), 'b--', 'linewidth', 2);                       % Nitrogen
        lh.lines(end+1) = semilogx(ax, m, rates.jLreal(ix), 'g--', 'linewidth', 2);                   % Light
        lh.lines(end+1) = semilogx(ax, m, rates.jDOC(ix), 'color', [181 100 30]/256, 'linewidth', 2); % DOC
        lh.lines(end+1) = semilogx(ax, m, rates.jMax(ix), 'k:');                                      % Max growth rate
        lh.lines(end+1) = semilogx(ax, m, rates.jTot(ix), 'k--', 'linewidth', 2);                     % Total growth
    end
end

hold(ax, 'off');

ylim(ax, [0 2]);
xlim(ax, calcXlim(p));
ylabel(ax, 'Gains (day^{-1})');
set(ax, 'XScale', 'log');

% Set legend captions depending on presence of silicate
if isfield(p, 'ixSi')
    cap = {'Feeding', 'N', 'Light', 'Si', 'DOC', 'Max. growth rate', 'Growth rate'};
else
    cap = {'Feeding', 'N', 'Light', 'DOC', 'Max. growth rate', 'Growth rate'};
end

lh.legend = legend(ax, cap, 'location', 'eastoutside', 'box', 'off');

end
