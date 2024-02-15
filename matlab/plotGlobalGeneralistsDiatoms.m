%
% Plot annual average biomass of generalists and diatoms
%
function plotGlobalGeneralistsDiatoms(sim, options)
arguments
    sim struct;
    options.sProjection string = 'fast';
end
%
% find generalists:
%
if sum(sim.p.typeGroups==1)
    typeGeneralist = 1;
else
    typeGeneralist = 5;
end
Bgeneralists = calcBiomassGroup(sim, typeGeneralist);
%
% Find diatoms:
%
if sum( sim.p.typeGroups==4 )
    typeDiatom = 4;
else
    typeDiatom = 3;
end
Bdiatoms = calcBiomassGroup(sim, typeDiatom);
%%
% plot
%
clf
tiledlayout(2,1)

% Total unicellular biomass:
nexttile
handle = panelGlobal(sim.x,sim.y,log10(Bgeneralists+Bdiatoms),[-2 2],sProjection=options.sProjection,...
    sUnits='g_C/m^2');

% diatom fraction:
nexttile
panelGlobal(sim.x, sim.y, Bdiatoms./(Bgeneralists+Bdiatoms),0:0.1:1,sProjection=options.sProjection,...
    sUnits='Fraction diatoms')


  
