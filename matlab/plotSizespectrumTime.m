%
% Plot a Sheldon size spectrum over time in a given water column
%
% In:
%  sim - simulation structure
%  iDepth - The index of the depth to use
%  lat, lon - latitude and longitude (only needed for a global simulation)
%
% Out:
%  s - structure with the Sheldon spectrum
%
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
        s.B = sim.B;
        sTitle = "Sheldon spectrum";

    case 'watercolumn'
        % Extract from a single water column:
        s.B = squeeze(sim.B(:,iDepth,:));
        sTitle = sprintf("Community sheldon spectrum at depth %3.0f m", sim.z(iDepth));

    case 'global'
        if isempty(lat)
            disp('Must specify latitude and longitude')
            stop
        else
            % Extract from global run:
            idx = calcGlobalWatercolumn(lat,lon,sim);
            s.B = squeeze(sim.B(:, idx.x, idx.y, iDepth, :));
            sTitle = sprintf("Sheldon spectrum at (%i,%i) and %3.0f m depth", [lat,lon,sim.z(iDepth)]);
        end

    otherwise
        disp('Wrong model type ',sim.p.nameModel)
        stop
end

% Create Sheldon size spectrum
for iTime = 1:length(sim.t)
    [mc, s.Bc(:,iTime)] = calcCommunitySpectrum(s.B, sim, iTime);
end

s.Bc(imag(s.Bc) ~=0) = 0.01;
s.Bc(s.Bc < 0.01) = 0.01;

clf
contourf( mc, sim.t, s.Bc', ...
    10.^linspace(-2,2,100),'linestyle','none')

axis tight
set(gca,'xscale','log')
set(gca, 'colorscale','log')

c = colorbar;
c.Label.String = "log_{10}({\mu}g_C/L)";
sgtitle(sTitle)
xlabel('Mass (\mugC)')
ylabel('Time (days)')

