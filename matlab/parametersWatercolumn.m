function p = parametersWatercolumn(p, nTMmodel)

arguments
    p struct
    nTMmodel {mustBeInteger} = 2; % Use 1.0 deg model as default
end

p = parametersGlobal(p,nTMmodel);

p.nameModel = 'watercolumn';

if nTMmodel==2
    p.bUse_parday_light=false;
end

p.tSave = 1;
