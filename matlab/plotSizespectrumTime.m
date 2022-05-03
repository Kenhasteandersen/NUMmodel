function s = plotSizespectrumTime(sim,iDepth,lat,lon)
% θ
arguments
    sim struct;
    iDepth {mustBeInteger} = 1;
    lat double = [];
    lon double = [];
end

m = sim.p.m(sim.p.idxB:end);

switch sim.p.nameModel
    
    case 'chemostat'
        s.B = sim.B;
        sTitle = "Sheldon spectrum";
    
    case 'watercolumn'
        % Extract from a single water column:
        s.B = squeeze(sim.B(iDepth,:,:))';
        sTitle = sprintf("Community sheldon spectrum at depth of max biomass: %3.0f m", sim.z(iDepth));
    
    case 'global'
    if isempty(lat)
        disp('Must specify latitude and longitude')
        stop
    else
        % Extract from global run:
        idx = calcGlobalWatercolumn(lat,lon,sim);
        s.B = squeeze(sim.B(idx.x, idx.y, iDepth, :, :))';
        sTitle = sprintf("Sheldon spectrum at %3.0f m", sim.z(iDepth));
    end

    otherwise
        disp('Wrong model type ',sim.p.nameModel)
        stop
end
    
% Create Sheldon size spectrum from biomasses by dividing with bin width:
for iTime = 1:length(sim.t)

    [mc, s.Bc(:,iTime)] = calcCommunitySpectrum(s.B, sim, iTime);

    % [mc, s.Bc(:,iTime)] = calcCommunitySpectrum(sim.p, s.B(:,iTime));
    %s.B(:,iTime) = s.B(:,iTime) ./ sim.p.mDelta(sim.p.idxB:end)' .* m';
end

s.Bc(imag(s.Bc) ~=0) = 0.01;
s.Bc(s.Bc < 0.01) = 0.01;

clf
% surface(mc, sim.t, log10(s.Bc)')

% s.Bc( s.Bc < 0.01 ) = 0.01;
% BB(imag(BB) ~= 0) = 0.0001;

contourf( mc, sim.t, s.Bc', ...
    10.^linspace(-2,2,100),'linestyle','none')

% surface(mc, sim.t, s.Bc')
% shading flat
axis tight
% caxis([-2,2])
set(gca,'xscale','log')
set(gca, 'colorscale','log')

c = colorbar;
c.Label.String = "log_{10}({\mu}g_C/L)";
title(sTitle)
xlabel('Mass (\mugC)')
ylabel('Time (days)')

