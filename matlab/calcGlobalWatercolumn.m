%
% Returns the indices of the watercolumn in a global simulation.
%
% In:
%  lat, lon - latitude and longitude
%  sim
%
% Out:
%  structure with fields 'x" and "y" and closest to lat and lon. Also includes the
%  depths in field "z".
%
function idx = calcGlobalWatercolumn(lat, lon, sim)

arguments
    lat, lon;
    sim struct;
end

if (lon<0)
    lon = lon+360;
end

dist=(sim.x*ones(1,length(sim.y))-lon).^2 + ((sim.y*ones(1,length(sim.x)))'-lat).^2;
shortest = min(dist(:));
ix = find(dist==shortest);
ix = ix(1);
idx.x = mod(ix, length(sim.x))+1;
idx.y = floor(ix/length(sim.x));
idx.z = find(sim.bathy(idx.x, idx.y,:)==1);
%
% Check whether we're asking for a land point:
%
if isempty(idx.z)
    error('The point %3.0f, %3.0f is on land.\n', lat,lon);
end

