function [m, Bc] = calcCommunitySpectrum(sim)

[m, idx] = sort(sim.p.m(3:end));
B = sim.B(end,idx);
Bcumm = cumsum(B);
Bc = gradient(Bcumm);