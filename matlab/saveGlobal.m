%
% Save the results of a simulation to be used for initial conditions
%
function saveGlobal(sim)

sim.N = sim.N(:,:,:,end);
sim.DOC = sim.DOC(:,:,:,end);
sim.B = sim.B(:,:,:,:,end);
sim.t = 0;
sim.p = sim.p;

save(sim.p.pathInit,'sim');

    
