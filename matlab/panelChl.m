%
% Makes a panel of the Chl in a water
function sim = panelChl(sim)

if sim.p.nameModel ~= 'watercolumn'
    error('Only implemented for watercolumn simulations\n');
end

if ~isfield(sim,'ChlVolume')
    sim = calcFunctions(sim);
end

surface(sim.t, -sim.z, log10(sim.ChlVolume'));
shading flat;
colorbar
xlabel('Time (days)')
ylabel('Depth (m)')
ylim([-200 0])

title('log_{10}(Chl) ({\mu}g_{Chl}/l)')