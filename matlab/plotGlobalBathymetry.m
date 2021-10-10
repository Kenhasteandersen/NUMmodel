%
% Plots the bathymetry
%
function plotGlobalBathymetry(sim, sProjection)

arguments
    sim struct;
    sProjection string = 'fast';
end

ixDepth = sum(sim.bathy,3);
depth = 0*ixDepth;
depth(ixDepth>0) = sim.z(ixDepth(ixDepth>0));

c = panelGlobal(sim.x,sim.y,depth,[0 5000],sTitle='Water depth',sProjection=sProjection);
c.Label.String  = 'm';

