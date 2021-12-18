%
% Plots the net primary production from a water column run
%
function func = plotWatercolumnFunctions(sim, options)

arguments
    sim struct;
    options.depthMax = 200;
end

for iTime = 1:length(sim.t)
    for k = 1:length(sim.z)
        if ~isnan(sim.N(k,iTime))
            % Get the functions per volume at each depth and time:
            u = [squeeze(sim.N(k,iTime)), ...
                squeeze(sim.DOC(k,iTime)), ...
                squeeze(sim.B(k,:,iTime))];
            [func.ProdGross(iTime,k), func.ProdNet(iTime,k),func.ProdHTL(iTime,k),func.eHTL(iTime,k),...
                func.Bpico(iTime,k),func.Bnano(iTime,k),func.Bmicro(iTime,k)] = ...
                getFunctions(u, sim.L(k,iTime), sim.T(k,iTime));
            
        end
    end
end

z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
dt = sim.t(2)-sim.t(1);
%t = [sim.t-0.5*dt sim.t(end)+0.5*dt];
%panelField(t, -z, func.ProdNet);
surface(sim.t,-z(1:end-1),func.ProdNet')
shading flat

xlabel('Time (days)')
ylabel('Depth (m)')
title('Net primary production (g_C/m^3/yr)')
axis tight
ylim([-options.depthMax, 0]);
caxis([0 20])
colorbar;

