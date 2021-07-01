%
% Plot a water column at latitude `lat`
% and longitude `lon` at time (days)
%
function plotGlobalWatercolumn(sim, time, lat,lon, bNewplot)

arguments
    sim struct;
    time double;
    lat double = [];
    lon double = [];
    bNewplot = true;
end

[~, iTime] = min(abs(sim.t-time));

if ~isempty(lat)
    % Extract water column from global simulation:
    idx = calcGlobalWatercolumn(lat,lon,sim);
    z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
    B = squeeze(double(sim.B(idx.x, idx.y, idx.z, :, iTime)));
    for i = 1:length(idx.z)
        u(i,:) = [sim.N(idx.x, idx.y, idx.z(i), iTime), ...
            sim.DOC(idx.x, idx.y, idx.z(i), iTime), B(i,:)];
        L(i) = sim.L(idx.x, idx.y, idx.z(i), iTime);
        T(i) = sim.T(idx.x, idx.y, idx.z(i), iTime);
    end
    
else
    % Extract data from water column simulation:
    z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
    B = squeeze(double(sim.B(:, :, iTime)));
    for i = 1:length(sim.z)
        u(i,:) = [sim.N(i, iTime), sim.DOC(i, iTime), B(i,:)];
        L(i) = sim.L(i,iTime);
        T(i) = sim.T(i,iTime);
    end
end
    
m = [sim.p.mLower(sim.p.idxB:end), sim.p.mLower(end)+sim.p.mDelta(end)];

if bNewplot
    clf
    tiledlayout(2,1,'tilespacing','compact','padding','compact')
end

%
% Biomass spectrum:
%
nexttile
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
for i = 1:length(z)-1
    if isfield(sim,'Si')
           rates = getRates(sim.p,[sim.N(idx.x, idx.y, idx.z(i), iTime), ...
            sim.Si(idx.x, idx.y, idx.z(i), iTime), ...
        sim.DOC(idx.x, idx.y, idx.z(i), iTime), B(i,:)],...
        sim.L(idx.x, idx.y, idx.z(i), iTime));
    else
        rates = getRates(sim.p, u(i,:), L(i), T(i));
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
    
