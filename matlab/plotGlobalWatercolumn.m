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
m = [sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];

clf
tiledlayout(2,1,'tilespacing','compact','padding','compact')
arguments
    sim struct;
    lat, lon double;
    iTime {mustBeInteger} = length(sim.t);
end

idx = calcGlobalWatercolumn(lat,lon,sim);
m = [sim.p.mLower(3:end), sim.p.mLower(end)+sim.p.mDelta(end)];
z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];

clf
tiledlayout(2,1,'tilespacing','compact','padding','compact')
%
% Biomass spectrum:
%
nexttile
B = squeeze(double(sim.B(idx.x, idx.y, idx.z, :, iTime)));
B(B<0) = 0;
panelField(m, -z, log10(B)');

set(gca,'xscale','log')
caxis([-5 2])

title(['Spectrum at lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')

cbar = colorbar;
cbar.Label.String  = 'log_{10} biomass (\mug C l^{-1})';
%
% Trophic strategy:
%
nexttile
for i = 1:length(idx.z)
    rates = getRates([sim.N(idx.x, idx.y, idx.z(i), iTime), sim.DOC(idx.x, idx.y, idx.z(i), iTime), B(i,:)],...
        sim.L(idx.x, idx.y, idx.z(i), iTime));
 %   [~, col] = calcTrophicStrategy(rates);
    for j=1:length(m)-1
        colStrategy(j,i,:) = ...
            [min(1, max(0, 3*(rates.jFreal(j)))), ...
            min(1, max(0, 3*(rates.jLreal(j)))), ...
            min(1, max(0, 3*(rates.jDOC(j))))];
  %      colStrategy(j,i,:) = col(j,:) * ;
    end
end

panelField(m,-z,colStrategy);
set(gca,'xscale','log')
xlabel('Mass (\mugC)')
ylabel('Depth (m)')
    
