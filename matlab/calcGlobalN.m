%
% Calculate the global N as a function of time.
% In:
%  sim struct from a global simulation
% Out:
%  N summed over all cells (in grammes).
%
function N = calcGlobalN(sim)

N = 0*sim.t;

load(sim.p.pathGrid, 'dv');

for i = 1:length(sim.t)
    tmp = sim.N(:,:,:,i).*dv/1000*1e-6;
    N(i) = sum(tmp(~isnan(tmp)));
    
    tmp = sim.B(:,:,:,:,i)/5.68;
    N(i) = N(i) + sum(tmp(~isnan(tmp)));
end
