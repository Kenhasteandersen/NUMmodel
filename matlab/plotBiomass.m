%
% Plot the biomass of each group (generalists, diatoms, active copepods and passive copepods)
%
% In:
%   sim - simulation strucutre (chemostat, watercolumn or global)
%   tDay - day to plot (end of the simulation by default)
%   Optinal:
%   options.depthMax - max depth for ylimit (for watercolumn and global)
%   options.pourcentage - true means the biomass values are in pourcetage (do the average of the % for each groups, use for plotComparedSimulation)
%                         false means the biomass values are in µgC/L (sum the biomass for each group)
%

function plotBiomass(sim,tDay,options)

arguments
    sim struct
    tDay = sim.t(end);
    options.depthMax = [];
    options.pourcentage = false;
end

[~, iTime] = min(abs(sim.t-tDay));
p = sim.p;
groups = p.typeGroups; groups(groups==100) = []; groups(groups==5) = 1; groups(groups==4) = 3;
[groups, iGroups] = unique(groups,'stable');
IX = cell(length(groups),2);

%
% Generalists
%
if ~isempty(find(groups==1,1))
    ixS=[p.ixStart(p.typeGroups==1) p.ixStart(p.typeGroups==5)];
    ixE=[p.ixEnd(p.typeGroups==1) p.ixEnd(p.typeGroups==5)];
    ix = [];
    for i = 1:length(ixS)
        ix = [ix ixS(i):ixE(i)];
    end
    IX{groups==1,1} = ix;
    IX{groups==1,2} = 'Generalists';
end
%
% Diatoms
%
if ~isempty(find(groups==3,1))
    ixS=[p.ixStart(p.typeGroups==3) p.ixStart(p.typeGroups==4)];
    ixE=[p.ixEnd(p.typeGroups==3) p.ixEnd(p.typeGroups==4)];
    ix = [];
    for i = 1:length(ixS)
        ix = [ix ixS(i):ixE(i)];
    end
    IX{groups==3,1} = ix;
    IX{groups==3,2} = 'Diatoms';
end
%
% Passiv Copepods
%
if ~isempty(find(p.typeGroups==10,1))
    ixS=p.ixStart(p.typeGroups==10);
    ixE=p.ixEnd(p.typeGroups==10);
    ix = [];
    for i = 1:length(ixS)
        ix = [ix ixS(i):ixE(i)];
    end
    IX{groups==10,1} = ix;
    IX{groups==10,2} = 'Passiv Copepods';
end
%
% Activ Copepods
%
if ~isempty(find(p.typeGroups==11,1))
    ixS=p.ixStart(p.typeGroups==11);
    ixE=p.ixEnd(p.typeGroups==11);
    ix = [];
    for i = 1:length(ixS)
        ix = [ix ixS(i):ixE(i)];
    end
    IX{groups==11,1} = ix;
    IX{groups==11,2} = 'Activ Copepods';
end

%
% Plotting
%
switch p.nameModel

    case 'chemostat'
        B=sim.B(iTime,:);
    case 'watercolumn'
        z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
        B=squeeze(sim.B(iTime,:,:));
        tiledlayout(1,length(groups),'tilespacing','compact','padding','compact')
    case 'global'
        tiledlayout(2,length(groups),'tilespacing','compact','padding','compact')
        Lat = [45.01 41.93 38.63 35.3 31.55 28.61 25.05 21.45 17.74 13.81 10.1 6.26 -1.15 -4.985 -8.765 -12.085 -15.035 -21.6 -24.905 -27.2 -30.01 -32.375 -34.58 -37.01 -39.21 -41.37 -43.79 -46.02];
        Lon = [-13.58 -16.02 -18.53 -20.92 -22.42 -24.36 -26.05 -27.75 -28.93 -28.29 -27.29 -26.38 -24.99 -24.97 -24.95 -24.93 -25.03 -25.05 -25.903 -28.25 -31.15 -33.67 -36.09 -38.83 -41.47 -44.01 -46.95 -49.87];
        Max = max(max(max(max(max(sim.B)))));
        Min = min(min(min(min(min(sim.B)))));
        mask=isnan(squeeze(sim.B(iTime,:,:,1,1)));

end
    
for i = 1:length(groups)

    ix = IX{i,1};
    [m, ord]=sort(p.m(ix));
    ix=ix-p.idxB+1;
    
    switch p.nameModel
        case 'chemostat'
            hold on
            plot(m,B(ix),'color',p.colGroup{iGroups(i)},'linewidth', 2.5)
            set(gca,'xscale','log')
            xlabel('Size (\mugC)')
            set(gca,'xtick',10.^(-9:2:5))
            ylabel('Biomass (\mugC/l)')
            hold off

        case 'watercolumn'
            nexttile
            contourf(m,-z,[B(1,ix(ord)); B(:,ix(ord))],'linestyle','none');
            set(gca,'xscale','log')%,'colorscale','log')
            xlabel('Size (\mugC)')
            set(gca,'xtick',10.^(-9:2:5))
            if groups(1)==groups(i)
                ylabel('Depth (m)')
            else
                set(gca,'yticklabel',[]);
            end
            if ~isempty(options.depthMax)
                ylim([-options.depthMax, 0]);
            end
            title(IX{i,2})
            cb = colorbar;

        case 'global'
            if options.pourcentage
                B = mean(sim.B(:,:,:,:,ix),5);
            else 
                B = squeeze(sum(sim.B(:,:,:,:,ix),5));
            end

            nexttile(length(groups)+i)
            [lat, lon] = panelGlobalTransect(sim,B,Lat,Lon,tDay,nbXticksMax=4,depthMax=options.depthMax);
            caxis([Min Max])
            colorbar off
            if groups(1)~=groups(i)
                set(gca,'yticklabel',[]);
            end

            nexttile(i)
            contourf(sim.x,sim.y,squeeze(B(iTime,:,:,1))','linestyle','none')
            xlabel('Latitude')
            if groups(1)==groups(i)
                ylabel('Longitude')
            else
                set(gca,'yticklabel',[]);
            end
            title(IX{i,2})
            caxis([Min Max])
            hold on
            contour(sim.x,sim.y,mask',[.1 .1],'k','LineWidth',0.1);
            plot(lon,lat,'k','linewidth',2)
            hold off
    end

end 

sgtitle(strcat('Biomass - day:',num2str(tDay)))

switch p.nameModel

    case 'chemostat'
        legend(IX{:,2})

    case 'global'
        cb = colorbar;
        cb.Label.String = 'Biomass ({\mug}C/l)'; 
        cb.Layout.Tile = 'east';

    case 'watercolumn'
        cb = colorbar;
        cb.Label.String = 'Biomass ({\mug}C/l)';
        cb.Layout.Tile = 'east';

end
