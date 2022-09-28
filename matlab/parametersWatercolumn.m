function p = parametersWatercolumn(p, nTMmodel)

arguments
    p struct
    nTMmodel {mustBeInteger} = 2; % Use hi-res ecco model as default
end

p = parametersGlobal(p,nTMmodel);

p.nameModel = 'watercolumn';

p.tSave = 1;
p.DiffBottom = 10; % Diffusivity of N out of the bottom m^2/day