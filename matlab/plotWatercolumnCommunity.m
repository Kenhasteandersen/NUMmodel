function plotWatercolumnCommunity(sim, time, lat,lon, options)

arguments
    sim struct;
    time = NaN;
    lat double = [];
    lon double = [];
    options.bNewplot  = true;
    options.depthMax {mustBePositive} = [];
end


if strcmp(sim.p.nameModel, 'global')


else
    z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
    % Extract data from water column simulation:
    if isnan(time)
        B = sim.B;
    else
        [~, iTime] = min(abs(sim.t-time));
        B = squeeze(double(sim.B(iTime,:, :))); % B(depth, grid)
    end

end

B(B<0) = 0;

if options.bNewplot
    clf
    tiles = tiledlayout(2,sim.p.nGroups,'tilespacing','compact','padding','compact');
    tiles.TileIndexing = 'columnmajor';
end

if isempty(options.depthMax)
    options.depthMax = max(z);
end

%
% Run over all depths
%

nDepth = size(B,2);

for iDepth = 1:nDepth

    if isnan(time)
        [mc, BB(iDepth,:)] = calcCommunitySpectrumWaterCol(squeeze(B(:,iDepth,:)), sim);
    else
        [mc, BB(iDepth,:)] = calcCommunitySpectrumWaterCol(squeeze(B(:,iDepth)), sim, iTime);
    end

end


BB = [BB(1,:); BB]; % Add dummy layer on top
BB( BB < 0.01 ) = 0.01;
% BB(imag(BB) ~= 0) = 0.0001;
%

if isnan(time)
    contourf( mc, -z, BB, ...
        10.^linspace(-2,0.5,100),'linestyle','none')
    caxis([0.01 4])
else
    contourf( mc, -z, BB, ...
        10.^linspace(-2,2,100),'linestyle','none')
    caxis([0.01 100])
end


% contourf( mc, -z, BB,'linestyle','none')

% set(gca,'xscale','log','colorscale','log')

xlabel('Mass (\mugC)')
ylabel('Depth (m)')

ylim([-200 0])

set(gca,'xscale','log','colorscale','log')
% set(gca,'xtick',10.^(-9:2), 'XTickLabel',[])

cbar = colorbar;
cbar.Label.String  = 'Biomass (\mug C l^{-1})';

% caxis([0.01 100])

if strcmp(sim.p.nameModel, 'watercolumn') && isnan(time)
    sgtitle(['Community Sheldon biomass ({\mu}gC/l) for the last ', num2str(sim.p.tEnd/2), ' days, ', 'lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176)])
elseif strcmp(sim.p.nameModel, 'watercolumn') && ~isnan(time)
    sgtitle(['Community Sheldon biomass ({\mu}gC/l) at day: ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176)])
end