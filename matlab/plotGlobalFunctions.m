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
tiles = tiledlayout(3,1,'TileSpacing','compact','padding','compact');

nexttile
cbar = panelGlobal(sim.x,sim.y,sim.ProdNetAnnual,'Net primary production (gC/m^2/yr)',sProjection);
cbar.Label.String = 'gC m{^-2} yr^{-1}';
%cbar.Visible='off';
%caxis([-3,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,sim.ProdHTLAnnual,'HTL production (gC/m^2/yr)',sProjection);
cbar.Label.String = 'gC m{^-2} yr^{-1}';
%cbar.Visible='off';
%caxis([-3,2])

nexttile
cbar = panelGlobal(sim.x,sim.y,sim.ProdHTLAnnual./sim.ProdNetAnnual,'eHTL',sProjection);
caxis([0,2])
%cbar.Location = 'SouthOutside';



