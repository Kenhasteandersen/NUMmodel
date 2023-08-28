% 
% Plot the compared simulation from compareSimulations
%
% In:
%   sim - the compared simulation (output of compareSimulations)
%   sim1 - one of the simulation used to get sim
%


function plotComparedSimulations(sim,sim1)

arguments
    sim struct;
    sim1 struct;
end


switch sim.p.nameModel

    case'chemostat'
        %
        % Nutrients
        %
        if isfield(sim,'Si')
            nutrients = {sim1.N, sim1.DOC,sim1.Si ; sim.N, sim.DOC, sim.Si};
        else
            nutrients = {sim1.N, sim1.DOC ; sim.N, sim.DOC};
        end

        figure
        tiledlayout(2,sim.p.nNutrients,'tilespacing','compact','padding','compact')
        for i = 1:length(nutrients)

            nexttile(i)
            plot(sim1.t,nutrients{1,i},'Color',sim1.p.colNutrients{i})
            title(sim.p.nameNutrientsLong{i})
            xticks([])

            nexttile(i+length(nutrients))
            plot(sim.t,nutrients{2,i},'Color',sim.p.colNutrients{i})
            title(strcat(sim.p.nameNutrientsShort{i},' variations with previous simulation (%)'))
            xlabel('time (day)')
        end

        
        %
        % Biomass
        %
        figure
        tiledlayout(2,1,'tilespacing','compact','padding','compact')
        nexttile
        plotBiomass(sim1);
        nexttile
        plotBiomass(sim);
        title('Varition with previous simulation (%)')
        
    case 'watercolumn'
        %
        % Localisation
        %
        if (sim.lon<0)
            sim.lon = sim.lon+360;
        end
        figure
        contourf(sim.x,sim.y,-squeeze(sim.bathy(:,:,1))','linestyle','none')
        colorbar('Ticks',[-1,0],...
         'TickLabels',{'Ocean','Land'})
        xlabel('longitude (°)')
        ylabel('latitude (°)')
        cmocean('haline',2,'pivot',0)
        hold on
        plot(sim.lon-3,sim.lat,'dy','markersize',7,'MarkerFaceColor','y')
        hold off
        text(sim.lon,sim.lat,strcat(num2str(sim.lat),'° N,',num2str(sim.lon),'° E'),'Color','y','fontsize',8)
        %
        % Biomass
        %
        figure
        plotBiomass(sim)
        sgtitle('Changes with previous simulation (%)')
        c=colorbar;
        c.Label.String ='Biomass Variations (%)';
        for i=1:4
            nexttile(i)
            clim([-100 100]);
            cmocean('balance','pivot',0)
        end
        %
        % Nutrients
        %
        if isfield(sim,'Si')
            nutrients = {sim1.N, sim1.DOC,sim1.Si ; sim.N, sim.DOC, sim.Si};
        else
            nutrients = {sim1.N, sim1.DOC ; sim.N, sim.DOC};
        end

        figure
        z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
        tiledlayout(2,sim.p.nNutrients,'tilespacing','compact','padding','compact')
        for i=1:length(nutrients)
            
            nexttile(i)
            plot([nutrients{1,i}(end,1) nutrients{1,i}(end,:)],-z,'Color',sim1.p.colNutrients{i})
            title(sim1.p.nameNutrientsLong{i})
            ylabel('Depth (m)')
            xlabel(strcat(sim.p.nameNutrientsShort{i},' (',sim.p.nameNutrientsUnits{i},')'))

            nexttile(i+length(nutrients))
            plot([nutrients{2,i}(end,1) nutrients{2,i}(end,:)],-z,'Color',sim.p.colNutrients{i})
            title(strcat(sim.p.nameNutrientsShort{i},' variations with previous simulation in %'))
            ylabel('Depth (m)')
            xlabel('Variation %')

        end

    case 'global'
        
        Lat = [45.01 41.93 38.63 35.3 31.55 28.61 25.05 21.45 17.74 13.81 10.1 6.26 -1.15 -4.985 -8.765 -12.085 -15.035 -21.6 -24.905 -27.2 -30.01 -32.375 -34.58 -37.01 -39.21 -41.37 -43.79 -46.02];
        Lon = [-13.58 -16.02 -18.53 -20.92 -22.42 -24.36 -26.05 -27.75 -28.93 -28.29 -27.29 -26.38 -24.99 -24.97 -24.95 -24.93 -25.03 -25.05 -25.903 -28.25 -31.15 -33.67 -36.09 -38.83 -41.47 -44.01 -46.95 -49.87];
        B = mean(sim.B,5);
        B1 = mean(sim1.B,5);
        nameLong = [sim.p.nameNutrientsLong(:)', {'Total Biomass'}];
        nameShort = [sim.p.nameNutrientsShort(:)', {'Biomass'}];
        units = [sim.p.nameNutrientsUnits(:)', {'{\mug}C/l'}];
        if isfield(sim,'Si')
            nutrients = {sim1.N, sim1.DOC,sim1.Si,B1 ; sim.N, sim.DOC,sim.Si,B};
        else
            nutrients = {sim1.N, sim1.DOC, B1 ; sim.N, sim.DOC, B};
        end
        s = {sim1,sim};
        
        for i=1:length(nutrients)
            
            figure
            tiledlayout(2,2,'tilespacing','compact','padding','compact')
            mask = isnan(squeeze(nutrients{1,i}(end,:,:,1)));
            s = [s;{strcat(nameShort{i},' ( ',units{i},')'),strcat(nameShort{i},' varation (%)')}];
            for j = 1:2
                
                nexttile(j+2)
                [lat,lon] = panelGlobalTransect(s{1,j},nutrients{j,i},Lat,Lon,s{1,j}.t(end),distMax=100,nbXticksMax=5);
                c=colorbar;
                c.Label.String =s{2,j};
                
                nexttile(j)
                contourf(s{1,j}.x,s{1,j}.y,squeeze(nutrients{j,i}(end,:,:,1))','linestyle','none')
                c=colorbar;
                c.Label.String =s{2,j};
                hold on
                contour(s{1,j}.x,s{1,j}.y,mask',[.1 .1],'k','LineWidth',0.1);
                plot(lon,lat,'k','linewidth',2)
                legend('','','Transect','Location','north')
                hold off

            end

            sgtitle(nameLong{i})
            
            nexttile(1); title('Latest simulation')
            
            nexttile(2); title('Variations with previous simulation (%)')
            clim([-100 100]);
            cmocean('balance','pivot',0)
            yticks([])
    
            nexttile(4)
            clim([-100 100]);
            cmocean('balance','pivot',0)
            yticks([])
            ylabel([])

            s(2,:) = [];
        end

        figure
        plotBiomass(sim,pourcentage=true)
        sgtitle('Variation of the Biomass with the previous Simulation (%)')
        groups = sim.p.typeGroups; groups(groups==100) = [];
        for i=1:length(unique(groups))*2
            nexttile(i)
            clim([-100 100]);
            cmocean('balance','pivot',0);
        end
        c = colorbar;
        c.Label.String = 'Biomass variation (%)';
        c.Layout.Tile = 'east';
end

if isfield(sim,'optionRates')
    plotComparedRates(sim,sim.rates)
end