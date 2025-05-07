function [s, tl] = plotSizeSpectrumTest(sim, time, iDepth, lat, lon, parent)
arguments
    sim struct;
    time = sim.p.tEnd;
    iDepth {mustBeInteger} = 1;
    lat double = [];
    lon double = [];
    parent {mustBeA(parent, {'matlab.ui.Figure', 'matlab.ui.container.Panel', 'matlab.ui.container.GridLayout', 'matlab.graphics.layout.TiledChartLayout'})} = []
end


% Create tiledlayout depending on parent
if isempty(parent)
    tl = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
elseif isa(parent, 'matlab.graphics.layout.TiledChartLayout')
    tl = parent; % Reuse the parent tiledlayout if already given
else
    tl = tiledlayout(parent, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
end

m = [sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
[~, iTime] = min(abs(sim.t - time));

switch sim.p.nameModel
    case 'chemostat'
        s.B = sim.B;
        if isfield(sim, 'Si')
            u = [sim.N(iTime), sim.DOC(iTime), sim.Si(iTime), squeeze(sim.B(iTime, :))];
        else
            u = [sim.N(iTime), sim.DOC(iTime), squeeze(sim.B(iTime, :))];
        end
        s.L = mean(sim.L);
        s.T = sim.T;

    case 'watercolumn'
        z = sim.z;
        s.B = squeeze(sim.B(:, iDepth, :));
        if isfield(sim, 'Si')
            u = [sim.N(iTime, iDepth), sim.DOC(iTime, iDepth), sim.Si(iTime, iDepth), ...
                squeeze(s.B(iTime, :))];
        else
            u = [sim.N(iTime, iDepth), sim.DOC(iTime, iDepth), squeeze(s.B(iTime, :))];
        end
        s.L = sim.L(iTime, iDepth);
        s.T = sim.T(iTime, iDepth);

    case 'global'
        idx = calcGlobalWatercolumn(lat, lon, sim);
        z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
        s.B = squeeze(sim.B(:, idx.x, idx.y, iDepth, :));
        if isfield(sim, 'Si')
            u = [sim.N(iTime, idx.x, idx.y, iDepth), ...
                sim.DOC(iTime, idx.x, idx.y, iDepth), ...
                sim.Si(iTime, idx.x, idx.y, iDepth), ...
                squeeze(sim.B(iTime, idx.x, idx.y, iDepth, :))'];
        else
            u = [sim.N(iTime, idx.x, idx.y, iDepth), ...
                sim.DOC(iTime, idx.x, idx.y, iDepth), ...
                squeeze(sim.B(iTime, idx.x, idx.y, iDepth, :))'];
        end
        s.L = sim.L(iTime, idx.x, idx.y, iDepth);
        s.T = sim.T(iTime, idx.x, idx.y, iDepth);
end

s.p = sim.p;
s.t = sim.t;

if ~isfield(sim, 'rates')
    sim.rates = getRates(sim.p, u, s.L, s.T);
end

s.rates = sim.rates;
sim_rates = sim.rates;
save('sim_rates.mat', 'sim_rates');

% --- Now plotting on the tiles ---
nexttile(tl)
panelSpectrum(s, iTime, bPlotStrategies=false);
xlabel('');
set(gca, 'XTickLabel', '');

nexttile(tl)
panelGains(sim.p, sim.rates);
set(gca, 'XTickLabel', '');

nexttile(tl)
panelLosses(sim.p, sim.rates);
set(gca, 'XTickLabel', '');
xlabel('');

nexttile(tl)
panelTrophicLevel(sim.p, s.B(iTime, :), sim.rates);

if strcmp(sim.p.nameModel, 'watercolumn')
    sgtitle(tl, ['Day = ', num2str(time), ...
                 ', lat = ', num2str(sim.lat), char(176), ...
                 ', lon = ', num2str(sim.lon), char(176), ...
                 ', depth of max biomass: ', num2str(z(iDepth)), ' m']);
end

end
