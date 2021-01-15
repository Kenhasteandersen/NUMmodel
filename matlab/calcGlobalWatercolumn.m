function idx = calcGlobalWatercolumn(lat, lon, sim)

if (lon<0)
    lon = lon+360;
end

dist=(sim.x*ones(1,length(sim.y))-lon).^2 + ((sim.y*ones(1,length(sim.x)))'-lat).^2;
shortest = min(dist(:));
ix = find(dist==shortest);
ix = ix(1);
idx.y = floor(ix/length(sim.x));
idx.x = mod(ix, length(sim.x));
idx.z = find( ~isnan([squeeze(sim.N(idx.x, idx.y,:,2))]) );

