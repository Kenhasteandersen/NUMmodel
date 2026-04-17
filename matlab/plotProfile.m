%
% Plots the profile of nutrients and all groups from a water column 
% simulation.
%
function plotProfile(sim, options)
arguments
    sim struct;
    %time double;
    %lat double = [];
    %lon double = [];
    options.bNewplot  = true;
    options.depthMax {mustBePositive} = 300;
end

p = sim.p;
z = sim.z;

if options.bNewplot
    clf
end

ixTime = floor(length(sim.t)/2):length(sim.t); % 

plot(mean(sim.N(ixTime,:),1), -z, linewidth=2, color=p.colNutrients{1})
hold on
plot(mean(sim.DOC(ixTime,:),1), -z,linewidth=2, color=p.colNutrients{2});
plot(mean(sim.Si(ixTime,:),1), -z,linewidth=2, color=p.colNutrients{3});

for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB+1;
    plot(mean(sum(sim.B(:,:,ix),3),1), -z, Color=p.colGroup{iGroup}, linewidth=2);
end

ylim([-options.depthMax,0])
set(gca,'xscale','log')
xlim([0.1,500])
legend([{'N'},{'DOC'},{'Si'},p.nameGroup], Location='southeast')
xlabel('Concentration ({\mu}g/l)')
ylabel('Depth (m)')