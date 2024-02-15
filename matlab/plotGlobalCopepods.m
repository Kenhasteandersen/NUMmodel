%
% Plot annual average biomass of active and passive copepods
%
function plotGlobalCopepods(sim, options)
arguments
    sim struct;
    options.sProjection string = 'fast';
end
p = sim.p;
%%
% Find active and passive copepods:
%
%Bcopepod = zeros(length(sim.t),128,64,15);
%Bgroup = zeros(2,128,64);
%typeGroup = 10:11;
%for j = 1:2
%    ixGroup = find(p.typeGroups==typeGroup(j));
%    for i = 1:length(ixGroup)
%        ix = (p.ixStart(ixGroup(i)):p.ixEnd(ixGroup(i))) - p.idxB+1;
%        Bcopepod = Bcopepod + sum( sim.B(:,:,:,:,ix),5 );
%    end
%    Bgroup(j,:,:) = calcIntegrateGlobal(sim,Bcopepod,true);
%end
%%
% Find geometric mean size:
%
% ixGroup = find(p.typeGroups==10 || p.typeGroups==11);
% ixStart = min( p.ixStart(ixGroup) ) - p.idxB+1;
% ixEnd = max ( p.ixEnd(ixGroup) ) - p.idxB+1;
% meanTime = mean(sim.B,1);
% for i = 1:length(sim.x)
%     for j = 1:length(sim.y)
%         for k = 1:length(sim.z)
            
Bactive = calcBiomassGroup(sim,11);
Bpassive = calcBiomassGroup(sim,10);

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


  
