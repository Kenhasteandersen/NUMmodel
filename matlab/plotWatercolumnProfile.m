
% Plots the profile of nutrients and all groups from a water column
% simulation.
%
% Averages over the last year.
%
function plotWatercolumnProfile(sim, options)
arguments
    sim struct;
    options.tDay = []; % Time of profile. Defaults to average over last year
    %time double;
    %lat double = [];
    %lon double = [];
    options.bNewplot  = true;
    options.depthMax {mustBePositive} = 300;
end

p = sim.p;
z = sim.z;

if isempty(options.tDay)
    %
    % Average over the last year:
    %
    ixStart = find(sim.t==sim.t(end)-365,1);
    if isempty(ixStart)
        ixStart=1;
    end
    ixTime = ixStart:length(sim.t); %
else
    ixTime = find(sim.t==options.tDay,1);
end
%
% Make the plot
%
if options.bNewplot
    clf
end

% Nutrients:
plot(mean(sim.N(ixTime,:),1), -z, linewidth=2, color=p.colNutrients{1})
hold on
plot(mean(sim.DOC(ixTime,:),1), -z,linewidth=2, color=p.colNutrients{2});
plot(mean(sim.Si(ixTime,:),1), -z,linewidth=2, color=p.colNutrients{3});

% Biomass groups:
for iGroup = 1:p.nGroups
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup)) - p.idxB+1;
    plot(mean(sum(sim.B(:,:,ix),3),1), -z, Color=p.colGroup{iGroup}, linewidth=1);
end

ylim([-options.depthMax,0])
set(gca,'xscale','log')
xlim([0.1,500])
legend([{'N'},{'DOC'},{'Si'},p.nameGroup], Location='southeast')
xlabel('Concentration ({\mu}g/l)')
ylabel('Depth (m)')
