%
% Import N and Si fields from WOA. Using data from 2018
%

addpath('~/Documents/Source/NUMmodel/matlab')
TMmodel = 1; % Which model to interpolate onto
%
% Use the annual average N field:
%
sFile = 'woa18_all_n00_01.nc'; % 

lat = ncread(sFile,'lat');
lon = ncread(sFile,'lon');
depth = ncread(sFile,'depth');

N = ncread(sFile,'n_an'); % Objectively analyzed mean fields for moles_concentration_of_nitrate_in_sea_water at standard depth levels
% micromoles_per_kilogram
molar_mass = 14.006747; % g/mol
N = N*molar_mass; % to micro g N / l
%
% Load our grid:
%
p = setupGeneralistsSimpleOnly; % dummy call just to get the grid
p = parametersGlobal(p,TMmodel);

sim = load(p.pathGrid,'x','y','z','dznom','bathy');

[X,Y,Z] = meshgrid(sim.x,sim.y,sim.z);
N = interp3(lat,lon+180,depth,N, Y,X,Z);
% Reorder dimensions:
N = rearrange(N,length(sim.x),length(sim.y),length(sim.z));
% Get rid of NANs:
N = replaceNAN(N);

save(p.pathN0,'N')
%
% Read the silicate file:
%
sFile = 'woa18_all_i00_01.nc';
molar_mass = 28.085;
Si = ncread(sFile,'i_an') * molar_mass; % Objectively analyzed mean fields for moles_concentration_of_nitrate_in_sea_water at standard depth levels
Si = interp3(lat,lon+180,depth,Si, Y,X,Z);
Si = rearrange(Si,length(sim.x),length(sim.y),length(sim.z));
Si = replaceNAN(Si);

save(p.pathSi0,'Si')

function A = rearrange(a,nx,ny,nz)
A = zeros(nx,ny,nz);
for i = 1:nx
    A(i,:,:) = squeeze(a(:,i,:));
end
end

function A = replaceNAN(a)
A = a;
% Replace NANs with means values at each depth level:
for j = 1:size(a,1)
    for i = 1:size(a,3)
        ixNAN = isnan(a(j,:,i));
        value = mean( a(j,~ixNAN,i) );
        if isnan(value)
            value = mean(a(~isnan(a(:))));
        end
        A(j,ixNAN,i) = value;
    end
end
end
