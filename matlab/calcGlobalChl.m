%
% Calculates the Chl from a global simulation in units of g Chl/m2
% See calcChl for algorithm
%
% In:
%   sim - simulation structure from simulateGlobal
%
% Out:
%   ChlArea - ChlA (g/m2)
%   ChlVolume - ChlA (g/volume)
%
function [ChlArea ChlVolume] = calcGlobalChl(sim)

ChlArea = zeros(length(sim.x), length(sim.y));
ChlVolume = zeros(length(sim.x), length(sim.y), length(sim.z));

for i = 1:length(sim.x)
    for j = 1:length(sim.y)
        for k = 1:length(sim.t)    
            for l = 1:length(sim.z)
                B = squeeze(sim.B(i,j,l,:,k))';
                if ~isnan(sum(B))
                    u = [squeeze(sim.N(i,j,l,k)), ...
                        squeeze(sim.DOC(i,j,l,k)), ...
                        B];
                    rates = getRates(sim.p, u, sim.L(i,j,l,k), sim.T(i,j,l,k));
                    tmp =  sum( calcChl( B, rates, sim.L(i,j,l,k))) /1000;
                    if ~isnan(tmp)
                        ChlArea(i,j) = ChlArea(i,j) + tmp * sim.dznom(l);
                        ChlVolume(i,j,l) = tmp;
                    end
                end
            end
        end
    end
end

ChlArea = ChlArea/length(sim.t);

