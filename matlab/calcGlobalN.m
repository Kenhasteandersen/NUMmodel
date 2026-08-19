%
% Calculate the global N as a function of time.
% In:
%  sim struct from a global simulation
% Out:
%  N summed over all cells (in grammes).
%
function N = calcGlobalN(sim)

S = inputRead;
rhoCN = S.general.rhoCN;

N = 0*sim.t;
N_BC = N;
p = sim.p;

load(sim.p.pathGrid, 'dv');

for i = 1:length(sim.t)
    % Nutrient in the N field:
    tmp = squeeze(sim.N(i,:,:,:)).*dv/1000*1e-6;
    N(i) = sum(tmp(~isnan(tmp)));
    % Nutrients in biomass gropus:
    for k = 1:p.n-p.nNutrients
        tmp = squeeze(sim.B(i,:,:,:,k)).*dv/1000*1e-6/rhoCN;
        N(i) = N(i) + sum(tmp(~isnan(tmp)));
    end
    % Nutrients diffusing up from the bottom BC:
    %N_BC(i) = p.dtTransport* ...
    %        p.BCmixing(p.idxN)./dzBottom'.*(BCvalue(:,p.idxN)-u(ixBottom, p.idxN));

    % Nutrients lost due to sinking POM:
    
end
