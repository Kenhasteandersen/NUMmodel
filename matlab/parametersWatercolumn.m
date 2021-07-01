function p = parametersWatercolumn(p, nTMmodel)

arguments
    p struct
    nTMmodel {mustBeInteger} = 2; % Use hi-res ecco model as default
end

p = parametersGlobal(p,nTMmodel);
p.tSave = 1;

