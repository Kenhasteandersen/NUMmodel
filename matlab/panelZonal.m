%
% Plots the zonal average (across longitudes) as a function of time
%
% field has dimensions time, x, y
%

function cbar = panelZonal(sim,field)

arguments
    sim struct;
    field (:,:,:);
end
%
% Average across latitudes
%
zonal = zeros(size(field,1), length(sim.y));
% Find non-land indexes:
ixDepth = sum(sim.bathy,3);

for i = 1:length(sim.y)
    % Average across non-land indexes:
    tmp = field(:,ixDepth(:,i)>0,i);
    
    zonal(:,i) = mean(tmp,2);
    
end
%
% Plot
%
surface(sim.t(end-size(field,1)+1:end), sim.y, zonal')
shading flat
set(gca,'colorscale','log')
clim([.01,1])
cbar = colorbar;
xlabel('Time (days)')
ylabel('Latitude')
axis tight