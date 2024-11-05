% 
% Plot the compared rates from compareRates (jLreal, jFreal & jDOC), for
% each group (generalists, diatoms, activ copepods and pasiv copepods)
%
% In:
%   sim - the compared simulation (output of compareSimulations)
%   sim1 - one of the simulation used to get sim
%

function plotComparedRates(sim,rates,options)

arguments
    sim struct
    rates struct
    options.depthMax = [];
end

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

J = {rates.jLreal,'jLreal'; rates.jFreal,'jFreal'; rates.jDOC,'jDOC'};

switch p.nameModel

    case 'chemostat'
        figure
        tiledlayout(length(J),1,'tilespacing','compact','padding','compact')
    case 'watercolumn'
        z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
        for k = 1:length(J)
            figure(500+k)
            tiledlayout(1,length(groups),'tilespacing','compact','padding','compact')
        end
    case 'global'
        for k = 1:length(J)
            figure(500+k)
            tiledlayout(2,length(groups),'tilespacing','compact','padding','compact')
        end
        Lat = [45.01 41.93 38.63 35.3 31.55 28.61 25.05 21.45 17.74 13.81 10.1 6.26 -1.15 -4.985 -8.765 -12.085 -15.035 -21.6 -24.905 -27.2 -30.01 -32.375 -34.58 -37.01 -39.21 -41.37 -43.79 -46.02];
        Lon = [-13.58 -16.02 -18.53 -20.92 -22.42 -24.36 -26.05 -27.75 -28.93 -28.29 -27.29 -26.38 -24.99 -24.97 -24.95 -24.93 -25.03 -25.05 -25.903 -28.25 -31.15 -33.67 -36.09 -38.83 -41.47 -44.01 -46.95 -49.87];
        mask=isnan(squeeze(sim.N(1,:,:,1)));

end

for i = 1:length(groups)

    ix = IX{i,1};
    [m, ord]=sort(p.m(ix));
    ix=ix-p.idxB+1;

    switch p.nameModel
        case 'chemostat'
            for k = 1:length(J)
                %j = J(k,:);
                nexttile(k)
                hold on
                plot(m,J{k,1}(ix),'color',p.colGroup{iGroups(i)},'linewidth', 1.5)
                set(gca,'xscale','log')
                xticks([])
                ylabel(J{k,2})
                hold off
            end

        case 'watercolumn'
            for k = 1:length(J)
                figure(500+k)
                nexttile
                contourf(m,-z,[J{k}(1,ix(ord)); J{k}(:,ix(ord))],'linestyle','none');
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
                colorbar;
                clim([-100 100]);
                cmocean('balance','pivot',0);
            end

        case 'global'

            for k = 1:length(J)

                j(1,:,:,:,:) = mean(J{k}(:,:,:,ix),4);
                
                figure(500+k)
                nexttile(length(groups)+i)
                [lat, lon] = panelGlobalTransect(sim,j,Lat,Lon,1,nbXticksMax=4,depthMax=options.depthMax);
                colorbar('off')
                clim([-100 100])
                cmocean('balance','pivot',0);
                if groups(1)~=groups(i)
                    set(gca,'yticklabel',[]);
                end
    
                nexttile(i)
                contourf(sim.x,sim.y,squeeze(j(1,:,:,1))',100,'linestyle','none')
                xlabel('Latitude')
                if groups(1)==groups(i)
                    ylabel('Longitude')
                else
                    set(gca,'yticklabel',[]);
                end
                title(IX{i,2})
                clim([-100 100])
                cmocean('balance','pivot',0);
                hold on
                contour(sim.x,sim.y,mask',[.1 .1],'k','LineWidth',0.1);
                plot(lon,lat,'k','linewidth',2)
                hold off
            end
    end

end 

%sgtitle(strcat('',num2str(tDay)))

switch p.nameModel

    case 'chemostat'
        lgd = legend(IX{:,2});
        xlabel('Size (\mugC)')
        set(gca,'xtick',10.^(-9:2:5))
        lgd.Layout.Tile = 'east';

    case 'global'
        for k = 1:length(J)
            figure(500+k)
            sgtitle(J{k,2})
            cb = colorbar;
            cb.Label.String = [J{k,2} ' variation (%)']; 
            cb.Layout.Tile = 'east';
        end

    case 'watercolumn'
        for k = 1:length(J)
            figure(500+k)
            sgtitle(J{k,2})
            cb = colorbar;
            cb.Label.String = [J{k,2} ' variation (%)'];
            cb.Layout.Tile = 'east';
        end

end