function s = plotSizespectrumTime(sim,iDepth,lat,lon)

arguments
    sim struct;
    iDepth {mustBeInteger} = 1;
    lat double = [];
    lon double = [];
end

m = sim.p.m(sim.p.idxB:end);

switch sim.p.nameModel
    
    case 'chemostat'
        s.B = sim.B';
    
    case 'watercolumn'
        % Extract from a single water column:
        s.B = squeeze(sim.B(iDepth,:,:));
    
    case 'global'
    if isempty(lat)
        disp('Must specify latitude and longitude')
        stop
    else
        % Extract from global run:
        idx = calcGlobalWatercolumn(lat,lon,sim);
        s.B = squeeze(sim.B(idx.x, idx.y, iDepth, :, :));
    end

    otherwise
        disp('Wrong model type ',sim.p.nameModel)
        stop
end
    
% Create Sheldon size spectrum from biomasses by dividing with bin width:
for iTime = 1:length(sim.t)
    [mc, s.B(:,iTime)] = calcCommunitySpectrum(m, s.B(:,iTime));
    %s.B(:,iTime) = s.B(:,iTime) ./ sim.p.mDelta(sim.p.idxB:end)' .* m';
end

clf
surface(mc, sim.t, log10(s.B)')
shading flat
axis tight
caxis([-5,2])
set(gca,'xscale','log')
colorbar
xlabel('mass')
ylabel('Time (days)')

