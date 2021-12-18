%
% Returns the indices of the watercolumn in a global simulation.
%
% In:
%  lat, lat - latitude and longitude
%
% Out:
%  structure with fields 'x" and "y" and closest to lat and lon. Also includes the
%  depths in field "z".
%
function idx = calcGlobalWatercolumn(lat, lon, sim)

arguments
    lat, lon {mustBeInteger};
    sim struct;
end

if (lon<0)
    lon = lon+360;
end

dist=(sim.x*ones(1,length(sim.y))-lon).^2 + ((sim.y*ones(1,length(sim.x)))'-lat).^2;
shortest = min(dist(:));
ix = find(dist==shortest);
ix = ix(1);
idx.x = mod(ix, length(sim.x));
idx.y = floor(ix/length(sim.x));
idx.z = find(sim.bathy(idx.x, idx.y,:)==1);

