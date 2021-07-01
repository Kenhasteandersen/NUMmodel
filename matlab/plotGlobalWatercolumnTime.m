%
% Plot a water column as a function of time at latitude `lat`
% and longitude `lon`. If lat and lon are empty it is assumed that the
% simulation is from "simulationWatercolumn".
%
function plotGlobalWatercolumnTime(sim, lat, lon, bNewPlot)

arguments
    sim struct;
    lat double = [];
    lon double = [];
    bNewPlot = true
end

if ~isempty(lat)
    % Extract the water column from a global run:
    idx = calcGlobalWatercolumn(lat,lon,sim);
    N = squeeze(double(sim.N(idx.x, idx.y, idx.z,:)));
    DOC = squeeze(double(sim.DOC(idx.x, idx.y, idx.z,:)));
    for i = 1:sim.p.nGroups
        B = squeeze((sum(double(sim.B(idx.x, idx.y, idx.z,...
            (sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1,:)),4)));
    end
    z = sim.z(idx.z);
else
    N = sim.N;
    DOC = sim.DOC;
    DOC(DOC<=0) = 1e-8;
    B = squeeze(sum(sim.B,2));
    z = sim.z + 0.5*sim.dznom;
end

t = sim.t;

if bNewPlot
    clf
    tiledlayout(2+isfield(sim,'Si')+sim.p.nGroups,1,'tilespacing','compact','padding','compact')
end

nexttile
surface(t,-z, N)
title(['Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%set(gca,'XTickLabel','');
%set(gca,'yscale','log')
set(gca,'ColorScale','log')
shading interp
axis tight
colorbar('ticks',10.^(-2:3))
caxis([0.1 100])

if isfield(sim,'Si')
    nexttile
    surface(t,-z, squeeze((double(sim.Si(idx.x, idx.y, idx.z,:)))))
    title(['Silicate, lat ', num2str(lat),', lon ', num2str(lon)])
    ylabel('Depth (m)')
    %set(gca,'XTickLabel','');
    %set(gca,'yscale','log')
    set(gca,'ColorScale','log')
    shading interp
    axis tight
    colorbar('ticks',10.^(-2:3))
    caxis([0.1 1000])
end

nexttile
surface(t,-z, DOC)
title(['DOC, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%xlabel('Concentration (\mugC l^{-1})')
set(gca,'colorscale','log')
%set(gca,'ColorScale','log')
shading interp
axis tight
colorbar
caxis([0.1,2])

for i = 1:sim.p.nGroups
    nexttile
    surface(t,-z, B)
    title(sim.p.nameGroup(i))
    ylabel('Depth (m)')
    xlabel('Time (days)')
    %set(gca,'yscale','log')
    set(gca,'ColorScale','log')
    shading interp
    axis tight
    colorbar
    caxis([0.1 100])
    colorbar('ticks',10.^(-2:3))
end

end