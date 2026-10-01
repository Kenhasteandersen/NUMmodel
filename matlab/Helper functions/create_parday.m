p = setupGeneralistsOnly(5);
p1 = parametersGlobal(p,1);
p2 = parametersGlobal(p,2);

load('~/Documents/NUMmodel/TMs/MITgcm_2.8deg/parday.mat');
sim1 = load(p1.pathGrid,'x','y','z','dznom','bathy'); % Get grid
sim2 = load(p2.pathGrid,'x','y','z','dznom','bathy'); % Get grid

[x,y] = meshgrid(sim1.y,sim1.x);
[xq,yq] = meshgrid(sim2.y, sim2.x);


for i = 1:730
    i
    p = squeeze(parday(:,:,1,i)); 
    ix = ~isnan(p) & p~=0;
    F = scatteredInterpolant(x(ix),y(ix),double(p(ix)),'nearest','nearest');

    parday2(:,:,1,i) = F(xq,yq);
end

parday = parday2;
save(p2.pathPARday,'parday')
%%
clf
tiledlayout(2,1)

nexttile
surface(squeeze(parday(:,:,1,1))'); shading flat

nexttile
surface(squeeze(parday2(:,:,1,1))'); shading flat