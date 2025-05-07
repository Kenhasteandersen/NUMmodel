function returnPlot= plotWatercolumnTimeTest(sim, lat, lon,plotType,appFigure, options)

arguments
    sim struct;
    lat double = [];
    lon double = [];
    plotType = [];
    appFigure  = []
    options.bNewPlot logical = true;
    options.depthMax {mustBePositive} = [];
    options.bOnlyLastYear = false;
    options.nLevels = 20;
end


switch sim.p.nameModel
    case 'global'
        idx = calcGlobalWatercolumn(lat,lon,sim);
        N = squeeze(double(sim.N(:,idx.x, idx.y, idx.z)))';
        DOC = squeeze(double(sim.DOC(:,idx.x, idx.y, idx.z)))';
        if isfield(sim,'Si')
            Si = squeeze(double(sim.Si(:,idx.x, idx.y, idx.z)))';
        end
        for i = 1:sim.p.nGroups
            B(i,:,:) = squeeze(sum(double(sim.B(:,idx.x, idx.y, idx.z, ...
                (sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1)),5))';
        end
        z = sim.z(idx.z);
    case 'watercolumn'
        N = sim.N';
        DOC = sim.DOC';
        if isfield(sim,'Si')
            Si = sim.Si';
        end
        for i = 1:sim.p.nGroups
            B(i,:,:) = squeeze(sum(sim.B(:,:,(sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1),3))';
        end
        z = sim.z + 0.5*sim.dznom;
        lat = sim.lat;
        lon = sim.lon;
    otherwise
        fprintf('Model type %s not supported.\n', sim.p.nameModel)
end

N(N<=0) = 1e-8;
DOC(DOC<=0) = 1e-8;

t = sim.t;
z = [0; z];
N = [N(1,:); N];
if isfield(sim,'Si')
    Si = max(0,[Si(1,:); Si]);
end
DOC = [DOC(1,:); DOC];
B(:,2:length(z),:) = B;
B(:,1,:) = B(:,2,:);

if isa(appFigure, 'matlab.graphics.layout.TiledChartLayout')
    % appFigure is already a tiledlayout
    disp('Using provided tiledlayout.');
    returnPlot = appFigure;
    delete(allchild(returnPlot));  % Clear previous plots
elseif isa(appFigure, {'matlab.ui.Figure', 'matlab.ui.container.Panel', 'matlab.ui.container.GridLayout'})
    % Create a new tiledlayout if given a Figure, Panel, or GridLayout
    if ~strcmp(plotType, 'Biomass')
        disp("No biomass seen, creating new 1x1 layout")
        returnPlot = tiledlayout(appFigure, 1, 1, 'TileSpacing', 'tight', 'Padding', 'loose');
    else
        disp("Biomass detected, creating new multi-row layout")
        returnPlot = tiledlayout(appFigure, sim.p.nGroups, 1, 'TileSpacing', 'tight', 'Padding', 'loose');
    end
else
    error('Invalid appFigure type provided.');
end


if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];
xlimit = [sim.t(1) sim.t(end)];
if options.bOnlyLastYear
    xlimit(1) = sim.t(end)-365;
end
switch plotType
    case 'Nitrogen'
        
        ax1 = nexttile(returnPlot);
        contourf(ax1, t, -z, N, logspace(-2, 3, options.nLevels), 'LineStyle', 'none');
        title(ax1, 'Nitrogen', 'FontWeight', 'normal');
        ylabel(ax1, 'Depth (m)');
        axis(ax1, 'tight');
        h = colorbar(ax1, 'ticks', 10.^(-2:3));
        h.Label.String = '{\mu}g_N/l';
        ylim(ax1, ylimit);
        xlim(ax1, xlimit);
        set(ax1, 'colorscale', 'log');
        set(ax1, 'XTickLabel', '');
        xlabel(ax1, 'Time (days)');

    case 'Dissolved silicate'
        if isfield(sim, 'Si')
            axSi = nexttile(returnPlot);
            contourf(axSi, t, -z, Si, logspace(-2, 3, options.nLevels), 'LineStyle', 'none');
            title(axSi, 'Silicate', 'FontWeight', 'normal');
            ylabel(axSi, 'Depth (m)');
            axis(axSi, 'tight');
            h = colorbar(axSi, 'ticks', 10.^(-2:3));
            h.Label.String = '{\mu}g_{Si}/l';
            ylim(axSi, ylimit);
            xlim(axSi, xlimit);
            set(axSi, 'colorscale', 'log');
            set(axSi, 'XTickLabel', '');
            xlabel(axSi, 'Time (days)');
        end

    case 'Dissolved organic matter'
        ax2 = nexttile(returnPlot);
        contourf(ax2, t, -z, DOC, logspace(-2, 2, options.nLevels), 'LineStyle', 'none');
        title(ax2, 'DOC', 'FontWeight', 'normal');
        ylabel(ax2, 'Depth (m)');
        axis(ax2, 'tight');
        h = colorbar(ax2, 'ticks', 10.^(-2:2));
        h.Label.String = '{\mu}g_C/l';
        ylim(ax2, ylimit);
        xlim(ax2, xlimit);
        set(ax2, 'colorscale', 'log');
        set(ax2, 'XTickLabel', '');
        xlabel(ax2, 'Time (days)');
        
    case 'Biomass'
        for i = 1:sim.p.nGroups
            nexttile
            %surface(t,-z, squeeze(B(i,:,:)))
            B(B < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
            contourf(t,-z,(squeeze(B(i,:,:))),[logspace(-2,3,options.nLevels)],'LineStyle','none')
            title( sim.p.nameGroup(i) ,'FontWeight','normal');
            ylabel('Depth (m)')
            axis tight
            h = colorbar('ticks',10.^(-2:3));
            h.Label.String = '{\mu}g_C/l';
            set(gca, 'colorscale','log')
            ylim(ylimit)
            xlim(xlimit)
            clim(10.^[-2,3])
            if i ~= sim.p.nGroups
                set(gca,'XTickLabel','')
            else
                xlabel('Time (days)')
            end
        end
    otherwise

        i = find(strcmp(sim.p.nameGroup, plotType));

        if isempty(i)
            error('Group name "%s" not found in sim.p.nameGroup.', plotType);
        end
        
        % Prepare the plot
        ax3 = nexttile(returnPlot);
        B(B < 0.01) = 0.01;  % Set small values to avoid plotting issues
        contourf(ax3, t, -z, squeeze(B(i, :, :)), logspace(-2, 3, options.nLevels), 'LineStyle', 'none');
        title(ax3, plotType, 'FontWeight', 'normal');
        ylabel(ax3, 'Depth (m)');
        axis(ax3, 'tight');

        h = colorbar(ax3, 'ticks', 10.^(-2:3));
        h.Label.String = '{\mu}g_C/l';

        set(ax3, 'colorscale', 'log');
        ylim(ax3, ylimit);
        xlim(ax3, xlimit);
        clim(ax3, [10^-2, 10^3]);
        xlabel(ax3, 'Time (days)');
        % if i ~= sim.p.nGroups
        %     set(ax3, 'XTickLabel', '');
        % else
        %     xlabel(ax3, 'Time (days)');
        % end




end



end

