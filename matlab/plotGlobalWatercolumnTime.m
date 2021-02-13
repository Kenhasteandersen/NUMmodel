%
% Plot a water column as a function of time at latitude `lat`
% and longitude `lon`.
%
function plotGlobalWatercolumnTime(lat,lon,sim)

idx = calcGlobalWatercolumn(lat,lon,sim);
t = sim.t;
z = sim.z(idx.z);

clf
subplot(3,1,1)
surface(t,-z, squeeze(log10(double(sim.N(idx.x, idx.y, idx.z,:)))))
title(['log10 Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%set(gca,'yscale','log')
shading interp
axis tight
colorbar

subplot(3,1,2)
surface(t,-z, squeeze(log10(double(sim.DOC(idx.x, idx.y, idx.z,:)))))
title(['log10 DOC, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%xlabel('Concentration (\mugC l^{-1})')
%set(gca,'yscale','log')
shading interp
axis tight
colorbar

subplot(3,1,3) 
surface(t,-z, squeeze(log10(sum(double(sim.B(idx.x, idx.y, idx.z,:,:)),4))))
title(['log10 Plankton, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
xlabel('Time (days)')
%set(gca,'yscale','log')
shading interp
axis tight
colorbar
caxis([-1 2])

end