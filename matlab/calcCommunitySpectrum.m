%
% Calculate the community size spectrum from all groups.
%
% In:
%  m - mass of all biomass groups
%  B - biomasses
%
% Out:
%  mc - masses of the community spectrum (sorted version of m)
%  Bc - Sheldon community spectrum
%
function [m, Bc] = calcCommunitySpectrum(m,B)

%[m, idx] = sort(sim.p.m(3:end));
%B = sim.B(end,idx);
[m, idx] = sort(m);
Bcumm = cumsum(B(idx));
Bc = gradient(Bcumm);