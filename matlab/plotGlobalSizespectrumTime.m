function s = plotGlobalSizespectrumTime(sim,iDepth,lat,lon)

arguments
    sim struct;
    iDepth {mustBeInteger} = 1;
    lat double = [];
    lon double = [];
end

m = sim.p.m(sim.p.idxB:end);

if ~isempty(lat)
    % Extract from global run:
    idx = calcGlobalWatercolumn(lat,lon,sim);
    s.B = squeeze(sim.B(idx.x, idx.y, iDepth, :, :));
else
    % Extract from a single water column:
    s.B = squeeze(sim.B(iDepth,:,:));
end

% Create Sheldon size spectrum from biomasses by dividing with bin width:
for iTime = 1:length(sim.t)
    s.B(:,iTime) = s.B(:,iTime) ./ sim.p.mDelta(sim.p.idxB:end)' .* m';
end

clf
surface(m, sim.t, log10(s.B)')
shading flat
axis tight
caxis([-5,2])
set(gca,'xscale','log')
colorbar
xlabel('mass')
ylabel('Time (days)')

