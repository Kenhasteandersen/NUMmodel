%
% Plot a water column at latitude `lat`
% and longitude `lon` at time (days)
%
function plotGlobalWatercolumn(sim, lat,lon, time, bNewplot)

arguments
    sim struct;
    lat, lon, time double;
    bNewplot = true;
end

idx = calcGlobalWatercolumn(lat,lon,sim);
m = [sim.p.mLower(sim.p.idxB:end), sim.p.mLower(end)+sim.p.mDelta(end)];
z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
[~, iTime] = min(abs(sim.t-time));

if bNewplot
    clf
    tiledlayout(2,1,'tilespacing','compact','padding','compact')
end

%
% Biomass spectrum:
%
nexttile
B = squeeze(double(sim.B(idx.x, idx.y, idx.z, :, iTime)));
B(B<0) = 0;
panelField(m, -z, (B)');

set(gca,'xscale','log','colorscale','log')


set(gca,'xtick',10.^(-9:2))
caxis([0.1 100])

title(['Spectrum at lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')

cbar = colorbar;
cbar.Label.String  = 'biomass (\mug C l^{-1})';
ylim([-200 0])
%
% Trophic strategy:
%
nexttile
for i = 1:length(idx.z)
    if isfield(sim,'Si')
            rates = getRates(sim.p,[sim.N(idx.x, idx.y, idx.z(i), iTime), ...
            sim.Si(idx.x, idx.y, idx.z(i), iTime), ...
        sim.DOC(idx.x, idx.y, idx.z(i), iTime), B(i,:)],...
        sim.L(idx.x, idx.y, idx.z(i), iTime));
    else
    rates = getRates(sim.p, [sim.N(idx.x, idx.y, idx.z(i), iTime), ...
        sim.DOC(idx.x, idx.y, idx.z(i), iTime), B(i,:)],...
        sim.L(idx.x, idx.y, idx.z(i), iTime));
    end
 %   [~, col] = calcTrophicStrategy(rates);
    for j=1:length(m)-1
        colStrategy(j,i,:) = ...
            [min(1, max(0, 6*(rates.jFreal(j)))), ...
            min(1, max(0, 3*(rates.jLreal(j)))), ...
            min(1, max(0, 3*(rates.jDOC(j))))];
  %      colStrategy(j,i,:) = col(j,:) * ;
    end
end

panelField(m,-z,colStrategy);
set(gca,'xscale','log')
set(gca,'xtick',10.^(-9:2))
xlabel('Cell mass (\mugC)')
ylabel('Depth (m)')
ylim([-200 0])
    
