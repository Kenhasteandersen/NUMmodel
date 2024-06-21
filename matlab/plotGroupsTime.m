%
% Plots the total biomass of each group from a chemostat simulation
%
function plotGroupsTime(sim)

p = sim.p;

if isnan(p.seasonalOptions.lat_lon) & p.seasonalOptions.seasonalAmplitude == 0
    sLegend = [];
    %
    % Biomass
    %
    for i = 1:p.nGroups
        ix = (p.ixStart(i):p.ixEnd(i)) - p.idxB + 1;
        B(i,:) = sum(sim.B( :, ix ),2);
        legendentries(i) = semilogy(sim.t, B(i,:), 'linewidth',2,'color',p.colGroup{i});
        sLegend{i} = p.nameGroup{i};
        hold on
    end
    %
    % Nutrients:
    %
    for i = 1:p.nNutrients
        legendentries(p.nGroups + i) = ...
            semilogy(sim.t, sim.u(:,i), 'color', p.colNutrients{i},'linewidth',2);
        sLegend{p.nGroups+i} = p.nameNutrientsLong{i};
    end

    hold off
    
    legend(legendentries, sLegend, 'location','northeastoutside','box','off')
    %ylim([0.1 max(B(:))])
    xlabel('Time (days)')
    ylabel('Total biomass (ugC/l)')

else
    %
    % Calculating the vectors of light and mixing rate for the right vector of
    % time (sim.t)
    %
    for i=1:size(sim.t,1)
        t_int = floor(mod(sim.t(i),365))+1;
        if t_int>365
            t_int = 365;
        end
        d(i) = p.d(t_int);
        L(i) = p.L(t_int);
    end
    
    %
    % Setup tiles:
    %
    clf
    tiledlayout(3,1,'tilespacing','compact','padding','compact')
    %
    % Biomass
    %
    nexttile
    for i = 1:p.nGroups
        ix = (p.ixStart(i):p.ixEnd(i)) - p.idxB + 1;
        B(i,:) = sum(sim.B( :, ix ),2);
        semilogy(sim.t, B(i,:), 'linewidth',2,'Color',p.colGroup{i});
        hold on
    end
    hold off
    legend(p.nameGroup,'box','off','location','southeast');
    set(gca,'XTickLabel','')
    ylabel('Total biomass (\mug_C/l)')
    xlim([min(sim.t) max(sim.t)])
    %
    % Mixing rate
    %
    nexttile
    semilogy(sim.t, d, 'linewidth', 2, 'color', 'r');
    set(gca,'XTickLabel','')
    ylabel('Mixing rate (1/day)')
    xlim([min(sim.t) max(sim.t)])
    %
    % Light
    %
    nexttile
    plot(sim.t, L, 'linewidth', 2, 'color', 'y');
    xlabel('Time (days)')
    ylabel('Light (\mu mol s^{-1} m^{-1})')
    xlim([min(sim.t) max(sim.t)])
end

