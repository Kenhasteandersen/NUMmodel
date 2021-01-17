function N = calcGlobalN(sim)

N = 0*sim.t;
for i = 1:length(sim.t)
    tmp = sim.N(:,:,:,i);
    N(i) = sum(tmp(~isnan(tmp)));
    
    tmp = sim.B(:,:,:,:,i)/5.68;
    N(i) = N(i) + sum(tmp(~isnan(tmp)));
end
