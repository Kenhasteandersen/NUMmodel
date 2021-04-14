%
% Plot a water column at latitude `lat`
% and longitude `lon` at time index iTime
%
function plotGlobalWatercolumn(sim, lat,lon, iTime)

arguments
    sim struct;
    lat, lon double;
    iTime {mustBeInteger} = length(sim.t);
end

idx = calcGlobalWatercolumn(lat,lon,sim);
z = sim.z(idx.z);

clf
B = squeeze(double(sim.B(idx.x, idx.y, :, :, iTime)));
ix = squeeze(~isnan( [sim.N(idx.x, idx.y, :, 12)] ));
imagesc( sim.p.mLower(3:end), sim.z(ix), log10(B(ix,:)) );
set(gca, 'xscale','log')
caxis([-5 2])

title(['Spectrum at lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
xlabel('Mass (\mugC)')

cbar = colorbar;
cbar.Label.String  = 'log_{10} biomass (\mug C l^{-1})';

axis tight
