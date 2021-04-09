function [tiles,sim] = plotGlobalFunctions(sim, sProjection)
arguments
    sim struct;
    sProjection string = 'fast';
end
%
% Get the functions:
%
sim = calcGlobalFunction(sim);

%%
clf
tiles = tiledlayout(3,1);%,'TileSpacing','compact');%,'padding','compact');
tiles.InnerPosition = [0.13,0.11,0.65,0.8150]; % Make space for colorbars

nexttile
cbar = panelGlobal(sim.x,sim.y,log10(sim.ProdNetAnnual(:,:,end)),'Net primary production',sProjection);
cbar.Label.String = 'log_{10}(gC m^{-2}yr^{-1})';
%cbar.Visible='off';
caxis([0,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,log10(sim.ProdHTLAnnual(:,:,end)),'HTL production',sProjection);
cbar.Label.String = 'log_{10}(gC m^{-2}yr^{-1})';
caxis([0,2])
%cbar.Visible='off';
%caxis([-3,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,sim.ProdHTLAnnual(:,:,end)./sim.ProdNetAnnual(:,:,end),'eHTL',sProjection);
caxis([0,1])
cbar.Label.String = '';
%cbar.Location = 'SouthOutside';


