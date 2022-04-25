function plotWatercolumnCommunity(sim, time, lat,lon, options)

arguments
    sim struct;
    time double;
    lat double = [];
    lon double = [];
    options.bNewplot  = true;
    options.depthMax {mustBePositive} = [];
end

[~, iTime] = min(abs(sim.t-time));

if strcmp(sim.p.nameModel, 'global')


else
    % Extract data from water column simulation:
    z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
    B = squeeze(double(sim.B(:, :, iTime))); % B(depth, grid)

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

nDepth = length(B(:,1));

for iDepth = 1:nDepth


    [mc, BB(iDepth,:)] = calcCommunitySpectrumWaterCol(B(iDepth,:), sim);

    set(gca,'xscale','log','colorscale','log')
    set(gca,'xtick',10.^(-9:2), 'XTickLabel',[])

    cbar = colorbar;
    cbar.Label.String  = 'Biomass (\mug C l^{-1})';
end


BB = [BB(1,:); BB]; % Add dummy layer on top
BB( BB<0.01 ) = 0.01;


contourf( mc, -z, BB, ...
    10.^linspace(-2,2,20),'linestyle','none')

set(gca,'xscale','log','colorscale','log')

xlabel('Mass (\mugC)')
ylabel('Depth (m)')

ylim([-200 0])


if strcmp(sim.p.nameModel, 'watercolumn')
    sgtitle(['Community Sheldon biomass ({\mu}gC/l) at day: ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176)])
end