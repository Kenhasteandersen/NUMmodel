%
% Plots the total biomass of each group from a chemostat simulation
%
function plotGroupsTime(sim)

p = sim.p;

for i = 1:p.nGroups
    ix = (p.ixStart(i):p.ixEnd(i)) - p.idxB + 1;
    B(i,:) = sum(sim.B( :, ix ),2);
    semilogy(sim.t, B(i,:), 'linewidth',2);
    hold on
end
hold off

legend(p.nameGroup);
ylim([0.1 max(B(:))])
xlabel('Time (days)')
ylabel('Total biomass (ugC/l)')

