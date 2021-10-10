%
% Plots the light at the surface (correction from the first cell).
%
function plotGlobalLight(sim, tDay, sProjection)

arguments
    sim struct;
    tDay double = sim.t(end);
    sProjection string = 'fast';
end

[~, iTime] = min(abs(sim.t-tDay));

L = squeeze(sim.L(:,:,1,iTime));
L = L*exp(sim.p.kw*sim.z(1)); % Correct from the first cell to the surface

c = panelGlobal(sim.x,sim.y,L,...
    sTitle='Daily average surface light', sProjection=sProjection);
c.Label.String  = '{\mu}E/m^2/s';
