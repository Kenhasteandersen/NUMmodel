%
% Plot the ecosystem functions from a global simulation. Calls
% "calcFunctions" to add the functions to the simulation structure. This
% call takes a long time, but the call is only made if the functions are
% not already present.
%
% In:
%  sim - simulation from a call to "simulateGlobal"
%
% Out:
%  sim - updated simulation strucure
%  tiles - the tiles with the plots
%
function [sim, tiles] = plotGlobalFunctions(sim, sProjection)
arguments
    sim struct;
    sProjection string = 'fast';
end
%
% Get the functions:
%
% sim = calcFunctions(sim);

%%
clf
set(gcf, 'Color', 'white')
tiles = tiledlayout(3,1);%,'TileSpacing','compact');%,'padding','compact');
tiles.InnerPosition = [0.13,0.11,0.65,0.8150]; % Make space for colorbars

nexttile
cbar = panelGlobal(sim.x,sim.y,log10(sim.ProdNetAnnual(end,:,:)),[1,3],...
    sTitle='Net primary production', sProjection=sProjection);
cbar.Label.String = 'log_{10}(mgC m^{-2}d^{-1})';
%cbar.Visible='off';
% caxis([1,3])
set(gca,'XTickLabel','')

nexttile
cbar = panelGlobal(sim.x,sim.y,log10(sim.ProdHTLAnnual(end,:,:)),[1,3],...
    sTitle='HTL production', sProjection=sProjection);
cbar.Label.String = 'log_{10}(mgC m^{-2}d^{-1})';
% caxis([1,3])

%cbar.Visible='off';
%caxis([-3,2])
set(gca,'XTickLabel','')

nexttile
cbar = panelGlobal(sim.x,sim.y,sim.ProdHTLAnnual(end,:,:)./sim.ProdNetAnnual(end,:,:),[0,1],...
    sTitle='\epsilon_{HTL}', sProjection=sProjection);
caxis([0,1])
cbar.Label.String = '';
set(gca,'XTickLabel','')

%cbar.Location = 'SouthOutside';



