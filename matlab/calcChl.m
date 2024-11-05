%
% Calculate ChlA in a size spectrum. Uses the relation from Edwards et al
% (2015) that gChl/gC are approximately proportional to the mass-specific light
% affinity.
%
% In:
%   simulation structure
%
% Out:
%   Total mass of ChlA in mg Chl/m3 and mg Chl/m2
%

% function [ChlVolume, ChlArea] = calcChl(sim)
%
% switch sim.p.nameModel
%
%     case 'chemostat'
%         ChlVolume = calculateChl(  sim.B(end,:), sim.rates.jLreal, sim.L );
%         ChlArea = ChlVolume * sim.p.depthProductiveLayer;
%
%     case 'global'
%         ChlArea = zeros(length(sim.x), length(sim.y));
%         ChlVolume = zeros(length(sim.x), length(sim.y), length(sim.z));
%
%         for i = 1:length(sim.x)
%             for j = 1:length(sim.y)
%                 for l = 1:length(sim.z)
%                     for k = 1:length(sim.t)
%                         B = squeeze(sim.B(i,j,l,:,k))';
%                         B(isnan(B)) = 0; % get rid of NaN
%                         if ~isnan(sim.N(i,j,l,k))
%                             u = [squeeze(sim.N(i,j,l,k)), ...
%                                 squeeze(sim.DOC(i,j,l,k)), ...
%                                 B];
%                             rates = getRates(sim.p, u, sim.L(i,j,l,k), sim.T(i,j,l,k));
%                             tmp =  calculateChl( B, rates, sim.L(i,j,l,k)) * 1000; % Convert to mg
%                             if ~isnan(tmp)
%                                 ChlArea(i,j) = ChlArea(i,j) + tmp * sim.dznom(l);
%                                 ChlVolume(i,j,l) = ChlVolume(i,j,l) + tmp;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         %
%         % Average over time:
%         %
%         ChlArea = ChlArea/length(sim.t);
%         ChlVolume = ChlVolume/length(sim.t);
% 
% end
% p.cA % correlation between (AL/m) vs. Chl-a:C in units mummol photons
% m^-2s^-1 day mugChl/mugC

function Bchl = calcChl(B,rates,L) %check units!
  cA = 1; % in units mummol photons m^-2s^-1 day mugChl/mugC
    if size(B)~=size(rates.jLreal)
%         error('Size of B is different from size of JLreal, should be the same');
        Bchl = cA*sum( B .* rates.jLreal' )/L; % in units of mu g Chl per l 
    else
        Bchl = cA*sum( B .* rates.jLreal )/L; % in units of mu g Chl per l 
    end
end 

%end
