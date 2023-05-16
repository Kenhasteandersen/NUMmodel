%
% Plot a global distribution of biomass and a water column in the same plot
%
% In:
%  time - time (in days)
%  lat, lon - latitude and longitude
%  Optional:
%  options.bNewplot - whether to clear the figure.
%  options.depthMax - max depth for ylimit.
%
function plotGlobalWatercolumn(sim, time, lat, lon, options)

arguments
    sim struct;
    time double;
    lat double = 60;
    lon double = -10;
    options.bNewplot  = true;
    options.depthMax {mustBePositive} = [];
end

[~, iTime] = min(abs(sim.t-time));

tiledlayout(4,1,'TileSpacing','compact','padding','compact');
%
% Global plot
%
nexttile(1,[2,1])

B = calcIntegrateGlobal(sim, sim.B(iTime,:,:,:,:));

c = panelGlobal(sim.x, sim.y, ...
    log10(B(:,:,1)), ...
    [0 2],...
    sProjection='eckert4');
c.Label.String  = 'g_C/m^2';

set(gca,'color', get(gcf,'color'));
gridm('off'); % Remove grid lines

plotm(lat,lon,'rp','markersize',20, 'MarkerFaceColor','red','markeredgecolor','red')

title('Total plankton biomass')
%annotation(figure1,'arrow',[0.5 0.5], [0.827571428571429 0.45]);

%
% Watercolumn:
%
nexttile
plotWatercolumn(sim, time, lat, lon, depthMax=200, bNewPlot=false);

nexttile(4)
%ylabel('')
%set(gca,'yticklabel','')
%xticklabel = get(gca,'xticklabel');
title('Trophic strategy')

nexttile(3)
%title('')
%set(gca,'xtick',10.^(-9:2),'xticklabel',xticklabel)
%xlabel('Cell mass (\mugC)')
%colorbar('off')
title('Size spectrum')

%annotation(gcf,'rectangle',...
%    [0.0277857142857143 0 0.9115 0.490476190476193],'edgecolor','red');

hold off

