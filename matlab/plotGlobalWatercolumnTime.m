%
% Plot a water column as a function of time at latitude `lat`
% and longitude `lon`.
%
function plotGlobalWatercolumnTime(lat,lon,sim, bNewPlot)

arguments
    lat double;
    lon double;
    sim struct;
    bNewPlot = true
end

idx = calcGlobalWatercolumn(lat,lon,sim);
t = sim.t;
z = sim.z(idx.z);


if bNewPlot
    clf
    tiledlayout(3,1,'tilespacing','compact','padding','compact')
end
nexttile
surface(t,-z, squeeze((double(sim.N(idx.x, idx.y, idx.z,:)))))
title(['Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%set(gca,'XTickLabel','');
%set(gca,'yscale','log')
set(gca,'ColorScale','log')
shading interp
axis tight
colorbar('ticks',10.^(-2:3))
caxis([0.1 100])

nexttile
surface(t,-z, squeeze((double(sim.DOC(idx.x, idx.y, idx.z,:)))))
title(['DOC, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%xlabel('Concentration (\mugC l^{-1})')
set(gca,'colorscale','log')
%set(gca,'ColorScale','log')
shading interp
axis tight
colorbar
caxis([0.1,2])

nexttile
surface(t,-z, squeeze((sum(double(sim.B(idx.x, idx.y, idx.z,:,:)),4))))
title(['Plankton, lat ', num2str(lat),', lon ', num2str(lon)])
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