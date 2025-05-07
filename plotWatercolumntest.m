function returnPlot = plotWatercolumn(sim, time, lat, lon, plotType, appFigure, options)

arguments
    sim struct;
    time double;
    lat double = [];
    lon double = [];
    plotType char = '';
    appFigure = [];
    options.bNewplot logical = true;
    options.depthMax {mustBePositive} = [];
end

[~, iTime] = min(abs(sim.t - time));

% Extract data depending on model type
switch sim.p.nameModel
    case 'global'
        idx = calcGlobalWatercolumn(lat, lon, sim);
        if isempty(idx.z)
            error('Not on land')
        end

        B = squeeze(double(sim.B(iTime, idx.x, idx.y, idx.z, :)));
        z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];

        for i = 1:length(idx.z)
            if isfield(sim, 'Si')
                u(i,:) = [sim.N(iTime, idx.x, idx.y, idx.z(i)), ...
                    sim.DOC(iTime, idx.x, idx.y, idx.z(i)), ...
                    sim.Si(iTime, idx.x, idx.y, idx.z(i)), ...
                    B(i,:)];
            else
                u(i,:) = [sim.N(iTime, idx.x, idx.y, idx.z(i)), ...
                    sim.DOC(iTime, idx.x, idx.y, idx.z(i)), ...
                    B(i,:)];
            end
            L(i) = sim.L(iTime, idx.x, idx.y, idx.z(i));
            T(i) = sim.T(iTime, idx.x, idx.y, idx.z(i));
        end

    case 'watercolumn'
        B = squeeze(double(sim.B(iTime, :, :)));
        z = [sim.z - 0.5*sim.dznom; sim.z(end) + 0.5*sim.dznom(end)];

        for i = 1:length(sim.z)
            if isfield(sim, 'Si')
                u(i,:) = [sim.N(iTime, i), sim.DOC(iTime,i), sim.Si(iTime,i), B(i,:)]; 
            else
                u(i,:) = [sim.N(iTime,i), sim.DOC(iTime,i), B(i,:)];
            end
            L(i) = sim.L(iTime, i);
            T(i) = sim.T(iTime, i);
        end
        lat = sim.lat;
        lon = sim.lon;

    otherwise
        error('Unsupported model type: %s', sim.p.nameModel)
end

B(B<0) = 0;

% Compute strategy and feeding level colors
nZ = length(z) - 1;
nGroups = sim.p.nGroups;
colStrategy = zeros(nGroups, nZ, 3);
colFeeding = colStrategy;
f = nan(nGroups, nZ);

for i = 1:nZ
    rates = getRates(sim.p, u(i,:), L(i), T(i));
    for j = 1:size(u,2) - sim.p.idxB + 1
        % Trophic strategy
        colStrategy(j,i,:) = [min(1, max(0, 6*rates.jFreal(j))), ...
            min(1, max(0, 3*rates.jLreal(j))), ...
            min(1, max(0, 3*rates.jDOC(j)))] ;

        % Feeding level
        f(j,i) = rates.jFreal(j)/rates.jFmaxx(j);
        if isnan(f(j,i))
            colFeeding(j,i,:) = [0 0 0];
        else
            fc = rates.jR/rates.jMax;
            if f(j,i) < fc
                colFeeding(j,i,:) = [0, 1, f(j,i)/fc];
            else
                colFeeding(j,i,:) = [min(1,3*f(j,i)), 0, 0];
            end
        end
    end
end

% Determine which groups to plot based on plotType
if isempty(plotType)
    groupIndices = 1:nGroups;  % Plot all groups
else
    groupIndices = find(strcmp(sim.p.nameGroup, plotType));  % Plot only selected group
end

% Adjust layout size dynamically based on the number of groups
nToPlot = numel(groupIndices);

% Set up the tiled layout (2, n) for multiple groups or (2, 1) for one group
if isa(appFigure, 'matlab.graphics.layout.TiledChartLayout') && isvalid(appFigure)
    returnPlot = appFigure;
    delete(allchild(returnPlot));  % Clear existing tiles only
else
    if isempty(appFigure)
        appFigure = figure;
    end
    returnPlot = tiledlayout(appFigure, 2, nToPlot, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
end

% Setup depth limits
if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];

% Plotting biomass spectra and strategy/feeding
Zmax = -Inf;
Zmin = Inf;

% Loop over the selected group indices
for k = 1:nToPlot
    iGroup = groupIndices(k);

    % Skip the group if it doesn't match plotType
    if ~isempty(plotType) && ~strcmp(sim.p.nameGroup{iGroup}, plotType)
        continue;
    end

    ix = sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup);
    m = [sim.p.mLower(ix), sim.p.mLower(ix(end))+sim.p.mDelta(ix(end))];
    ixB = ix - sim.p.idxB + 1;
    Delta = sim.p.mUpper(ix) ./ sim.p.mLower(ix);

    % Normalized biomass
    BB = B(:, ixB) ./ log(Delta);
    BB = [BB(1, :); BB];
    BB(BB < 0.01) = 0.01;

    if size(BB, 2) == 1
        mm = [sim.p.mLower(ix(1)), sim.p.mUpper(ix(1))];
        BB(:, 2) = BB(:, 1);
    else
        mm = sim.p.m(ix);
    end

    % Top panel: biomass
    ax1 = nexttile(returnPlot, k);  % Assign the correct tile
    contourf(ax1, mm, -z, BB, 10.^linspace(-2, 2, 20), 'linestyle', 'none');
    xlim(ax1, [min(mm), max(mm)]);
    set(ax1, 'xscale', 'log', 'colorscale', 'log');
    set(ax1, 'xtick', 10.^(-9:2), 'XTickLabel', []);
    title(ax1, sim.p.nameGroup{iGroup}, 'FontWeight', 'normal');
    if k == 1
        ylabel(ax1, 'Depth (m)');
    else
        set(ax1, 'yticklabel', []);
    end
    ylim(ax1, ylimit);
    caxis(ax1, [0.01 100]);

    Zmax = max(Zmax, max(BB(:)));
    Zmin = min(Zmin, min(BB(:)));

    if k == nToPlot
        cbar = colorbar(ax1);
        cbar.Label.String = 'Sheldon biomass ({\mu}g C/l)';
        if Zmin == Zmax
            Zmax = Zmin + 0.0001;
        end
        caxis(ax1, [Zmin Zmax]);
    end

    % Bottom panel: trophic strategy / feeding
    ax2 = nexttile(returnPlot, k + nToPlot);  % Assign the bottom tile
    xlim(ax2, [min(mm), max(mm)]);
    if sim.p.typeGroups(iGroup) < 10
        panelField(m, -z, colStrategy(ixB, :, :), ax2);
    else
        panelField(m, -z, colFeeding(ixB, :, :), ax2);
    end
    set(ax2, 'xscale', 'log');
    if sim.p.typeGroups(iGroup) < 10 || sim.p.typeGroups(iGroup) >= 100
        set(ax2, 'xtick', 10.^(-9:2:5));
        xlabel(ax2, 'Cell size ({\mu}gC)');
    else
        set(ax2, 'xtick', 10.^(-9:1:5));
        xlabel(ax2, 'Body size ({\mu}gC)');
    end
    if k == 1
        ylabel(ax2, 'Depth (m)');
    else
        set(ax2, 'yticklabel', []);
    end
    ylim(ax2, ylimit);
end
