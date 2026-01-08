%
% Calculate ChlA in a size spectrum. Uses the relation from Edwards et al
% (2015) that gChl/gC are approximately proportional to the mass-specific light
% affinity.
%
% In:
%   simulation structure
%
% Out:
%   Total mass of ChlA in mug Chl/liter
%

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
