%
% Integrate a global field over depth. If the field is 5-dimensional, it is
% assumed to be a biomass spectrum, that is then summed.
%
% In:
%  sim: simulation structure
%  field: field to integrate, e.g., field = sim.B or field = sim.N
%  bAverageTime: whether to average over time (default false)
%
% Out:
%  Three dimensional field in units of g/m2
%
function field = calcIntegrateGlobal(sim, field, bAverageTime)

arguments
    sim struct;
    field double;
    bAverageTime = false;
end

field(isnan(field)) = 0;

if length(size(field)==5) % Assume a spectrum
    field = sum(field,5); % Total biomass in the group
end
% Integrate over depth:
dz = sim.dznom;
field = squeeze( sum(field.*reshape(dz ,1,1,1,numel(dz)),4) / 1000); % g/m2

if bAverageTime
    field = mean(field,1);
end