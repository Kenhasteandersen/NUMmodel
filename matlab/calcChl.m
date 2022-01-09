%
% Calculate ChlA in a size spectrum. Uses the relation from Edwards et al
% (2015) that gChl/gC are approximately proportional to the mass-specific light
% affinity. I calculate the effective affinity as jLreal/L * 3, where the 
% factor 3 is a fit-by-eye.
%
% In:
%   B - biomasses of size groups
%   rates - a "rates" structure from getRates
%   L - Light (in mu mol photons per m2 per s)
%
% Out:
%   Total mass of ChlA (mg Chl/m2)
%
function BChl = calcChl(sim)

switch sim.p.nameModel
    
    case 'chemostat'
        BChl = 3 * sum( sim.B(end,:) .* sim.rates.jLreal' )/sim.L ...
             *sim.p.depthProductiveLayer; % Unit conversion
end

