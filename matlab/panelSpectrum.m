function lh = panelSpectrum(sim, ixTime, ax, options)
arguments
    sim struct;
    ixTime {mustBeInteger} = length(sim.t); % Defaults to last time step
    ax = [];
    options.bPlotStrategies = true;
end

p = sim.p;

% Activate axis
if isempty(ax)
    ax = gca;
end
axes(ax); % Explicitly use given axis

cla(ax); % Clear current content
hold(ax, 'on');

if options.bPlotStrategies
    rates = sim.rates;
    % Background color depending on trophic strategies
    [strategy, col] = calcTrophicStrategy(p, rates,ax);
end

sLegend = {};

% Plot spectra background
for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix - p.idxB + 1;

    % Averaging over last half of simulation
    ixAve = find(sim.t > sim.t(end)/2);

    if length(sim.t) > 1
        Blower = min(sim.B(ixAve, ixB)) ./ log(Delta);
        Bupper = max(sim.B(ixAve, ixB)) ./ log(Delta);

        patch(ax, ...
            [m, m(end:-1:1)], [Blower, Bupper(end:-1:1)], ...
            p.colGroup{iGroup}, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end
end

set(ax, 'XScale', 'log', 'YScale', 'log');

% Community spectrum
[mc, Bc] = calcCommunitySpectrum(sim.B, sim);
legendentries(1) = loglog(ax, mc, Bc, 'LineWidth', 4.5, 'Color', [0.7, 0.7, 0.7]);
sLegend{1} = 'Community spectrum';

% Group spectra
for iGroup = 1:p.nGroups
    ix = p.ixStart(iGroup):p.ixEnd(iGroup);
    m = p.m(ix);
    Delta = p.mUpper(ix)./p.mLower(ix);
    ixB = ix - p.idxB + 1;

    % Avoid log(0)
    sim.B(sim.B <= 0) = 1e-100;

    % Average spectrum
    B = exp(mean(log(sim.B(ixAve, ixB) ./ log(Delta)), 1));

    legendentries(iGroup+1) = loglog(ax, m, B, 'LineWidth', 2, 'Color', p.colGroup{iGroup});

    % Spectrum at specific time step
    B = sim.B(ixTime, ixB) ./ log(Delta);
    loglog(ax, m, B, ':', 'LineWidth', 1, 'Color', p.colGroup{iGroup});

    sLegend{iGroup+1} = p.nameGroup{iGroup};
end

ylim(ax, [0.0001, 500]);
xlim(ax, calcXlim(sim.p));

xlabel(ax, 'Mass ({\mu}g_C)');
ylabel(ax, 'Sheldon biomass ({\mu}g_C/L)');

lh = legend(ax, legendentries, sLegend, 'Box', 'off');
lh.Location = 'eastoutside';
end
