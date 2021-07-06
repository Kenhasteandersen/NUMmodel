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
%   Total mass of ChlA in each size group (mug Chl/l)
%
function BChl = calcChl(B, rates, L)

BChl = 3 * B .* rates.jLreal'/L;

