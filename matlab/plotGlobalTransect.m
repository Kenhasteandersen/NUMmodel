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
%   options.bAMTtrack - plots the AMT track (does not require that Lat and
%                       Lon are set)
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
    options.bAMTtrack = true;
end

if options.bAMTtrack
    Lat = [45.01 41.93 38.63 35.3 31.55 28.61 25.05 21.45 17.74 13.81 10.1 6.26 -1.15 -4.985 -8.765 -12.085 -15.035 -21.6 -24.905 -27.2 -30.01 -32.375 -34.58 -37.01 -39.21 -41.37 -43.79 -46.02];
    Lon = [-13.58 -16.02 -18.53 -20.92 -22.42 -24.36 -26.05 -27.75 -28.93 -28.29 -27.29 -26.38 -24.99 -24.97 -24.95 -24.93 -25.03 -25.05 -25.903 -28.25 -31.15 -33.67 -36.09 -38.83 -41.47 -44.01 -46.95 -49.87];
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
panelGlobalTransect(sim,sim.N,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax);
c=colorbar;
c.Label.String =strcat(sim.p.nameNutrientsShort(1),' ( ',sim.p.nameNutrientsUnits(1),')');
xlabel('')
set(gca,'XTickLabel','');
%
%DOC
%
nexttile
panelGlobalTransect(sim,sim.DOC,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax);
c=colorbar;
c.Label.String =strcat(sim.p.nameNutrientsShort(2),' ( ',sim.p.nameNutrientsUnits(2),')');
xlabel('')
set(gca,'XTickLabel','');
%
%Silicium
%
if isfield(sim,'Si')
    nexttile
    panelGlobalTransect(sim,sim.Si,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax);
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
panelGlobalTransect(sim,B,Lat,Lon,tDay,distMax=options.distMax,depthMax=options.depthMax);
c=colorbar;
c.Label.String ='Biomass ({\mu}C/l';

if options.bAMTtrack
    sgtitle('Approximate AMT track - average over 1 year')
end
