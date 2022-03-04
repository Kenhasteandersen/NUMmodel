%
% Plots a profile of the watercolumn at "day"
% 
function plotWatercolumnProfile(sim, day)

arguments
    sim struct;
    day double = sim.t(end);
end

clf
%
% Extract fields:
%
iTime = find(sim.t >= day, 1);
N = sim.N(:,iTime);
L = sim.L(:,iTime);
DOC = sim.DOC(:,iTime);
B = sum(sim.B(:,:,iTime),2);
T = sim.T(:,iTime);
z = sim.z;
%
% Make panel:
%
plot(N, -z,'b','linewidth',2)
hold on
plot(L, -z, 'y','linewidth',2)
plot(DOC,-z,'color',[165 42 42]/256,'linewidth',2)
plot(B,-z,'linewidth',2)
plot(T,-z,'r','linewidth',2)
set(gca,'xscale','log')

ylabel('Depth (m)')
xlim([0.01 1000])

legend({'N ()','Light ()','DOC','Biomass ({\mu}g_C/l)','T'}, ...
    'location','southoutside')