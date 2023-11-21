%
% Calculate the global N as a function of time.
% In:
%  sim struct from a global simulation
% Out:
%  N summed over all cells (in grammes).
%
function N = calcGlobalN(sim)

S = inputRead;
rhoCN = S.input_general.rhoCN;

N = 0*sim.t;

load(sim.p.pathGrid, 'dv');

for i = 1:length(sim.t)
    tmp = squeeze(sim.N(i,:,:,:)).*dv/1000*1e-6;
    N(i) = sum(tmp(~isnan(tmp)));
    
    for k = 1:sim.p.n-sim.p.nNutrients
        tmp = squeeze(sim.B(i,:,:,:,k)).*dv/1000*1e-6/rhoCN;
        N(i) = N(i) + sum(tmp(~isnan(tmp)));
    end
    
end
