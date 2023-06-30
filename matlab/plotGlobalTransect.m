%
% Plot nutrients and biomass along transect defined by Lat and Lon (using pannelGlobalTransect)
%
% In:
%   sim - Global simulation structure
%   tDay - time (in days)
%   Lat,Lon - latitude and longitude
%
% If tDay<0, plots the average over the last year of the simulation
%

function []=plotGlobalTransect(sim,tDay,Lat,Lon)

arguments
    sim struct;
    tDay double;
    Lat double = [];
    Lon double = [];
end

clf
if isfield(sim,'Si')
    tiledlayout(4,1,'tilespacing','compact','padding','compact')
else
    tiledlayout(3,1,'tilespacing','compact','padding','compact')
end
%
%Nitrogen
%
nexttile
panelGlobalTransect(sim,tDay,sim.N,Lat,Lon,distMax=100)
c=colorbar;
c.Label.String =strcat(sim.p.nameNutrientsShort(1),' ( ',sim.p.nameNutrientsUnits(1),')');
xlabel('')
set(gca,'XTickLabel','');
%
%DOC
%
nexttile
panelGlobalTransect(sim,tDay,sim.DOC,Lat,Lon,distMax=100,depthMax=500)
c=colorbar;
c.Label.String =strcat(sim.p.nameNutrientsShort(2),' ( ',sim.p.nameNutrientsUnits(2),')');
xlabel('')
set(gca,'XTickLabel','');
%
%Silicium
%
if isfield(sim,'Si')
    nexttile
    panelGlobalTransect(sim,tDay,sim.Si,Lat,Lon,distMax=100,depthMax=2000)
    c=colorbar;
    c.Label.String =strcat(sim.p.nameNutrientsShort(3),' ( ',sim.p.nameNutrientsUnits(3),')');
    xlabel('')
    set(gca,'XTickLabel','');
end
%
%Biomass
%
B=sum(sim.B,5);
nexttile
panelGlobalTransect(sim,tDay,B,Lat,Lon,distMax=100,depthMax=500)
c=colorbar;
c.Label.String ='Biomass ({\mu}C/l';
