%
% Plot a transect from the global simulation
%
% In:
%   sim - simulation structure
%   field - data to plot, must be a 4-dimensional matrix with (time, longitude, latitude, depth)
%   Lat,Lon - latitude and longitude
%   tDay - time (in days), equals to -1 by default
%          If tDay<0, plots the average over the last year of the simulation
%   Optional:
%   options.depthMax - max depth for ylimit
%   options.distMax - approximate distance between two plotting points (in km), 100 km by default
%                     (adds points beetewen the points defined by Lat and Lon)
%   options.nbXticksMax - nomber of ticks display on the x axes (longitude), 10 by default
%
% Out:
%   lat,lon - exact coordinates of the plotted transect
%  
%
% To change the color to log scale do : set(gca,'colorscale','log')
%
%
function [lat,lon]=panelGlobalTransect(sim,field,Lat,Lon,tDay,options)


arguments
    sim struct;
    field (:,:,:,:);
    Lat double = [];
    Lon double = [];
    tDay double = -1;
    options.depthMax {mustBePositive} = [];
    options.distMax {mustBePositive} = 100;
    options.nbXticksMax {mustBePositive} = 10;
end

lat=Lat;
lon=Lon;
%
% Adds Points
%
if ~isempty(options.distMax)
    dist=distance(Lat(1:end-1),Lon(1:end-1),Lat(2:end),Lon(2:end)).*pi/180*6371;
    nPoints=round(dist/options.distMax);
    nPoints(nPoints==0)=1;
    nPoints=nPoints-1;
    for i=find(nPoints>0)
        x=(1:nPoints(i)).*((Lon(i+1)-Lon(i))/(nPoints(i)+1))+Lon(i);
        y=(1:nPoints(i)).*((Lat(i+1)-Lat(i))/(nPoints(i)+1))+Lat(i);
        lon=[lon(1:i+sum(nPoints(1:i-1))) x lon(sum(nPoints(1:i-1))+i+1:end)];
        lat=[lat(1:i+sum(nPoints(1:i-1))) y lat(sum(nPoints(1:i-1))+i+1:end)];
    end
end
%
% 
%
xyz=table('Size',[length(lat) 3],'VariableTypes',["double" "double" "cell"],'VariableNames',["x" "y" "z"]);
for i=1:length(lat)
    idx=calcGlobalWatercolumn(lat(i), lon(i), sim);
    xyz(i,:)=struct2table(idx,"AsArray",true);
end
idx=unique(xyz(:,1:2),"stable");
idx=table2struct(idx,"ToScalar",true);
idx.z=1:length(sim.z);
clear xyz;


%
% Extract data from field:
%
z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
data=zeros(length(idx.x),length(idx.z));
coor=cell(length(idx.x),1);

for i=1:length(idx.x)
    if tDay<0
        %Average over the last year
        if sim.t(end)<365 %if the simulation lasted less than one year then the average is done over all the simulation
            data(i,:)=mean(squeeze(field(:,idx.x(i),idx.y(i),idx.z)));
        else
            deltat=sim.t(2)-sim.t(1);
            idxT=round(365/deltat);
            data(i,:)=mean(squeeze(field(end-idxT+1:end,idx.x(i),idx.y(i),idx.z)));
        end
    else
        [~, iTime] = min(abs(sim.t-tDay));
        data(i,:)=squeeze(field(iTime,idx.x(i),idx.y(i),idx.z));
    end
    coor{i}=[num2str(round(sim.y(idx.y(i)),1)) '°N, ' num2str(round(sim.x(idx.x(i)),1)) '°E'];
end

data=[data(:,1) data]; % Add dummy layer on top
%
%Plotting
%
contourf(1:length(idx.x),-z(:,1),data','linestyle','none')
colorbar;
xlabel('Coordinates')
ylabel('Depth (m)')

if length(idx.x) > options.nbXticksMax
    index = round(2:length(idx.x)/options.nbXticksMax:length(idx.x));
    xticks(index) 
    xticklabels(coor(index))
else
    xticks(1:length(idx.x))
    xticklabels(coor)
end

if ~isempty(options.depthMax)
    ylim([-options.depthMax, 0]);
end

%
% Mask
%
mask=isnan(data);
hold on
contour(1:length(idx.x),-z(:,1),mask',[.01 .01],'k','LineWidth',0.1)
hold off 

lon = sim.x(idx.x);
lat = sim.y(idx.y);