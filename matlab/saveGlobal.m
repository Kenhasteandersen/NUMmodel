%
% Save the results of a simulation to be used for initial conditions
%
function saveGlobal(sim, path)

arguments
    sim struct
    path string = sim.p.pathInit;
end

savesim.N = sim.N(end,:,:,:);
savesim.DOC = sim.DOC(end,:,:,:);
if isfield(sim, 'Si')
    savesim.Si = sim.Si(end,:,:,:);
end
sim.B(sim.B<0) = 0;
savesim.B = sim.B(end,:,:,:,:);
savesim.t = 0;
savesim.p = sim.p;
siminit = savesim;

save(path,'siminit','-v7.3');

    
