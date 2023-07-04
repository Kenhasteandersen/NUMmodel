%
% Plot nutrients and biomass along transect defined by Lat and Lon (using pannelGlobalTransect)
%
% In:
%   sim - Global simulation structure
%   Lat,Lon - latitude and longitude
%   tDay - time (in days), equals to -1 by default
%          If tDay<0, plots the average over the last year of the simulation
% Optional:
%   options.depthMax - max depth for ylimit
%   options.distMax - approximate distance between two plotting points (in km)
%                     (adds points beetewen the points defined by Lat and Lon)
%
%
function []=plotGlobalTransect(sim,Lat,Lon,tDay,options)

arguments
    sim struct;
    Lat double = [];
    Lon double = [];
    tDay double = -1;
    options.depthMax {mustBePositive} = [];
    options.distMax {mustBePositive} = 100;
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
panelGlobalTransect(sim,sim.N,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax)
c=colorbar;
c.Label.String =strcat(sim.p.nameNutrientsShort(1),' ( ',sim.p.nameNutrientsUnits(1),')');
xlabel('')
set(gca,'XTickLabel','');
%
%DOC
%
nexttile
panelGlobalTransect(sim,sim.DOC,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax)
c=colorbar;
c.Label.String =strcat(sim.p.nameNutrientsShort(2),' ( ',sim.p.nameNutrientsUnits(2),')');
xlabel('')
set(gca,'XTickLabel','');
%
%Silicium
%
if isfield(sim,'Si')
    nexttile
    panelGlobalTransect(sim,sim.Si,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax)
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
panelGlobalTransect(sim,B,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax)
c=colorbar;
c.Label.String ='Biomass ({\mu}C/l';
