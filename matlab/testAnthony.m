function returnGrid = testAnthony(sim, lat, lon, appContainer, options)
% testAnthony creates UI axes in a grid inside the provided UI container.
% appContainer: a uipanel or UIFigure content area where plots are placed.
%
% The options structure and simulation processing remains unchanged.

arguments
    sim struct;
    lat double = [];
    lon double = [];
    appContainer=[];  % Pas de type ici
    options.bNewPlot logical = true;
    options.depthMax {mustBePositive} = [];
    options.bOnlyLastYear logical = false;
    options.nLevels double = 20;
end

% Vérification manuelle
if ~isa(appContainer, 'matlab.ui.container.Panel') && ~isa(appContainer, 'matlab.ui.Figure')
    error('L''argument appContainer doit être un Panel ou une Figure.');
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

if options.bNewPlot
    clf
    returnPlot = tiledlayout(appContainer,2+isfield(sim,'Si')+sim.p.nGroups,1, ...
        'tilespacing','tight','padding','tight');
    
end

if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];
xlimit = [sim.t(1) sim.t(end)];
if options.bOnlyLastYear
    xlimit(1) = sim.t(end)-365;
end

%% Nitrogen
ax1 = nexttile(returnPlot);
contourf(ax1,t,-z,N,logspace(-2,3,options.nLevels),'LineStyle','none')
title(ax1,'Nitrogen','FontWeight','normal')
ylabel(ax1, 'Depth (m)')
axis(ax1, 'tight')
h = colorbar(ax1,'ticks',10.^(-2:3));
h.Label.String = '{\mu}g_N/l';
ylim(ax1, ylimit)
xlim(ax1, xlimit)
set(ax1, 'colorscale','log')
set(ax1, 'XTickLabel','')
ax1.ButtonDownFcn = @(src, event) onTileClicked(ax1, event, 'Nitrogen',sim);
%% Silicate (si présent)
if isfield(sim,'Si')
    axSi = nexttile(returnPlot);
    contourf(axSi, t, -z, Si, logspace(-2,3,options.nLevels), 'LineStyle', 'none')
    title(axSi, 'Silicate','FontWeight','normal')
    ylabel(axSi, 'Depth (m)')
    axis(axSi, 'tight')
    h = colorbar(axSi, 'ticks', 10.^(-2:3));
    h.Label.String = '{\mu}g_{Si}/l';
    ylim(axSi, ylimit)
    xlim(axSi, xlimit)
    set(axSi, 'colorscale', 'log')
    set(axSi, 'XTickLabel', '')
end
    axSi.ButtonDownFcn = @(src, event) onTileClicked(src, event, 'Silicate',sim);
%% DOC
ax2 = nexttile(returnPlot);
contourf(ax2,t,-z,DOC,logspace(-2,2,options.nLevels),'LineStyle','none')
title(ax2,'DOC','FontWeight','normal')
ylabel(ax2, 'Depth (m)')
axis(ax2, 'tight')
h = colorbar(ax2,'ticks',10.^(-2:2));
h.Label.String = '{\mu}g_C/l';
ylim(ax2, ylimit)
xlim(ax2, xlimit)
set(ax2, 'colorscale','log')
set(ax2, 'XTickLabel','')
    ax2.ButtonDownFcn = @(src, event) onTileClicked(src, event, 'DOC',sim);
%% Biomasse pour chaque groupe
for i = 1:sim.p.nGroups
    ax3 = nexttile(returnPlot);
    B(B < 0.01) = 0.01;
    contourf(ax3,t,-z,squeeze(B(i,:,:)),[logspace(-2,3,options.nLevels)],'LineStyle','none')
    title(ax3, sim.p.nameGroup(i) ,'FontWeight','normal');
    ylabel(ax3, 'Depth (m)')
    axis(ax3, 'tight')
    h = colorbar(ax3,'ticks',10.^(-2:3));
    h.Label.String = '{\mu}g_C/l';
    set(ax3, 'colorscale','log')
    ylim(ax3, ylimit)
    xlim(ax3, xlimit)
    clim(ax3, [10.^-2, 10.^3])

    if i ~= sim.p.nGroups
        set(ax3, 'XTickLabel','')
    else
        xlabel(ax3, 'Time (days)')
    end
    ax3.ButtonDownFcn = @(src, event) onTileClicked(src, event, sim.p.nameGroup(i),sim);
end

%% Titre général
if strcmp(sim.p.nameModel, 'watercolumn')
    sgtitle(appContainer, ['Water column at lat = ', num2str(lat), char(176), ', lon = ', num2str(lon), char(176)])
end

end
function onTileClicked(src, ~, titleText,sim,options)
% Get the X coordinate of the click
close all;
clickX = floor(src.CurrentPoint(1,1)); % Extract X coordinate from click

% Display the X coordinate in the command window
disp(['Tile clicked: ', titleText, ', X = ', floor(num2str(clickX))]);

figure
plotWatercolumn(sim,clickX,'depthMax',200);
% Change the background color of the clicked axes to highlight it
src.Color = [0.9, 0.9, 0.9]; % Light gray

figure
% Find the depth of maximum biomass:
Bdepth = sum( sim.B(clickX,:,:),3 );
iDepth = find(Bdepth==max(Bdepth));

%plotSizespectrum(sim,options.tDayPlot,iDepth);
plotSizeSpectrumTest(sim,clickX,iDepth);
% plotSizespectrum(sim,iDepth);

figure
% Find the depth of maximum average biomass:
Bdepth = sum(sum( sim.B(:,:,:),3 ),1);
iDepth = find(Bdepth==max(Bdepth));
plotSizespectrumTime(sim,iDepth);
end
