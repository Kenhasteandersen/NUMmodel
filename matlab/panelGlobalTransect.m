%
% Plot a transect from the global simulation
%
% In:
%   sim - simulation structure
%   tDay - time (in days)
%   field - data to plot
%   Lat,Lon - latitude and longitude
%   Optional:
%   options.depthMax - max depth for ylimit
%   options.distMax - approximativ distance max between two plotting points (in km)
%
% To change the color to log scale do : set(gca,'colorscale','log')
%
% If tDay<0, plots the average over the last year of the simulation
%
function []=panelGlobalTransect(sim,tDay,field,Lat,Lon,options)


arguments
    sim struct;
    tDay double;
    field (:,:,:,:);
    Lat double = [];
    Lon double = [];
    options.depthMax {mustBePositive} = [];
    options.distMax {mustBePositive} = [];
end

lat=Lat;
lon=Lon;
%
%Adds Points
%
if ~isempty(options.distMax)
    dist=distance(Lat(1:end-1),Lon(1:end-1),Lat(2:end),Lon(2:end)).*pi/180*6371;
    nPoints=round(dist/options.distMax);
    nPoints(nPoints==0)=1;
    nPoints=nPoints-1;
    for i=find(nPoints>0)
        x=(1:nPoints(i)-1)*((Lon(i+1)-Lon(i))/nPoints(i))+Lon(i);
        y=(1:nPoints(i)-1)*((Lat(i+1)-Lat(i))/nPoints(i))+Lat(i);
        lon=[lon(1:i+sum(nPoints(1:i-1))) x lon(sum(nPoints(1:i-1))+1:end)];
        lat=[lat(1:i+sum(nPoints(1:i-1))) y lat(sum(nPoints(1:i-1))+1:end)];
    end
end
%
%
%
xyz=table('Size',[length(lat) 3],'VariableTypes',["double" "double" "cell"],'VariableNames',["x" "y" "z"]);
for i=1:length(lat)
    idx=calcGlobalWatercolumn(lat(i), lon(i), sim);
    idx.z(sim.bathy(idx.x, idx.y,:)==0)=find((sim.bathy(idx.x, idx.y,:)==0));
    xyz(i,:)=struct2table(idx,"AsArray",true);
end
idx=unique(xyz(:,1:2),"stable");
idx=table2struct(idx,"ToScalar",true);
idx.z=cell2mat(table2array(xyz(1,3)));
clear xyz;
%
% Extract data from field:
%
z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
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
xticks(1:length(idx.x))
xticklabels(coor)
xlabel('Coordinates')
ylabel('Depth (m)')
if ~isempty(options.depthMax)
    ylim([-options.depthMax, 0]);
end
