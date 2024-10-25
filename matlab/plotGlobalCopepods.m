%
% Plot annual average biomass of active and passive copepods
%
function plotGlobalCopepods(sim, options)
arguments
    sim struct;
    options.sProjection string = 'fast';
end
            
Bactive = calcBiomassGroup(sim,12);
Bpassive = calcBiomassGroup(sim,13);

%%
% plot
%
clf
tiledlayout(2,1)

% Total copepod biomass:
nexttile
handle = panelGlobal(sim.x,sim.y,log10(Bactive+Bpassive),[-2 2],sProjection=options.sProjection,...
    sUnits='g_C/m^2');


% Active/passive:
nexttile
panelGlobal(sim.x, sim.y, Bactive./(Bactive+Bpassive),[0 1],sProjection=options.sProjection,...
    sUnits='Fraction active')


  
