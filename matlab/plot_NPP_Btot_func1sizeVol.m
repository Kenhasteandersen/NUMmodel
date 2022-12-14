% plots annual average NPP (gC/m2/yr)
% Total biomass (mugC/L) at surface
% Annual avg Chl-a(gCm-3)
% for sim produced by calcFunctions1sizeVol.m
sim5000calF=simGDfun;
%%
figure(55)
clf(55)
subplot(1,3,1)
surface(sim5000calF.ProdNetAnnual'); shading flat; colorbar; title("annual average NPP (gC/m2/yr)")
caxis([0 5000]); axis tight

subplot(1,3,2)
Btotal=squeeze(sum(sim5000calF.B(:,:,1,:,:),4)); %sum over size classes at surface
Btotal_av=mean(Btotal,3);
surface(log10(Btotal_av)'); shading flat; colorbar; title("Total biomass (\mugC/L) at surface")
caxis([0 3]);
axis tight


%%
% for Chl I will only take the mean of the last 12 months at surface
% figure(7)
% clf(7)
Chl_a=squeeze(mean(sim5000calF.ChlVolume(:,:,:,1),1));
subplot(1,3,3)
surface(log10(Chl_a)'); shading flat; colorbar; title("Annual average Chl-a (gCm^{-3})")
% caxis([-2.5 1]);
axis tight


