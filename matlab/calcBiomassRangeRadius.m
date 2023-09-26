%
% Return the biomass within the size range r1:r2.
% Also output the fraction of biomass from each size
% group that falls into the range.
%
function [Brange] = calcBiomassRangeRadius(B,r, rmin, rmax)
    ix = find(r>=rmin & r<rmax);
    Brange = sum(B(ix));