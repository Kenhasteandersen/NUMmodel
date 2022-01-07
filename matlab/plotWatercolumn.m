%
% Plot a water column from either a global or a water column simulation.
%
% Warning: the "setup" function needs to be called before this plot is
% made. If not matlab may crash or the results be incorrect.
%
% In:
%  time - time (in days)
%  lat, lon - latitude and longitude (only for global simulation)
%  Optional:
%  options.bNewplot - whether to clear the figure.
%  options.depthMax - mx depth for ylimit.
%
function plotWatercolumn(sim, time, lat,lon, options)

arguments
    sim struct;
    time double;
    lat double = [];
    lon double = [];
    options.bNewplot  = true;
    options.depthMax {mustBePositive} = [];
end

[~, iTime] = min(abs(sim.t-time));

if strcmp(sim.p.nameModel, 'global')
    % Extract water column from global simulation:
    idx = calcGlobalWatercolumn(lat,lon,sim);
    if isempty(idx.z)
        error('Not on land')
    end
    z = [sim.z(idx.z)-0.5*sim.dznom(idx.z); sim.z(idx.z(end))+0.5*sim.dznom(idx.z(end))];
    B = squeeze(double(sim.B(idx.x, idx.y, idx.z, :, iTime)));
    for i = 1:length(idx.z)
        if isfield(sim,'Si')
            u(i,:) = [sim.N(idx.x, idx.y, idx.z(i), iTime), ...
                sim.DOC(idx.x, idx.y, idx.z(i), iTime), ...
                sim.Si(idx.x, idx.y, idx.z(i), iTime), ...
                B(i,:)];
        else
            u(i,:) = [sim.N(idx.x, idx.y, idx.z(i), iTime), ...
                sim.DOC(idx.x, idx.y, idx.z(i), iTime), B(i,:)];
        end
        L(i) = sim.L(idx.x, idx.y, idx.z(i), iTime);
        T(i) = sim.T(idx.x, idx.y, idx.z(i), iTime);
    end
    
else
    % Extract data from water column simulation:
    z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
    B = squeeze(double(sim.B(:, :, iTime)));
    for i = 1:length(sim.z)
        if isfield(sim,'Si')
            u(i,:) = [sim.N(i, iTime), sim.DOC(i, iTime), ...
                sim.Si(i,iTime), B(i,:)];
        else
            u(i,:) = [sim.N(i, iTime), sim.DOC(i, iTime), B(i,:)];
        end
        L(i) = sim.L(i,iTime);
        T(i) = sim.T(i,iTime);
    end
end

B(B<0) = 0;
%
% Calclulate trophic strategy:
%
colStrategy = zeros(length(u)-sim.p.idxB+1, length(z)-1, 3);
for i = 1:length(z)-1
    rates = getRates(sim.p, u(i,:), L(i), T(i));
    for j=1:size(u,2)-sim.p.idxB+1
        colStrategy(j,i,:) = ...
            [min(1, max(0, 6*(rates.jFreal(j)))), ...
             min(1, max(0, 3*(rates.jLreal(j)))), ...
             min(1, max(0, 3*(rates.jDOC(j))))];
    end
end
%
% Setup tiles:
%
if options.bNewplot
    clf
    tiles = tiledlayout(2,sim.p.nGroups,'tilespacing','compact','padding','compact');
    tiles.TileIndexing = 'columnmajor';
    nexttile
end

if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];
%
% Run over all groups.
%
for iGroup = 1:sim.p.nGroups
    ix = sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup);
    m = [sim.p.mLower(ix), sim.p.mLower(ix(end))+sim.p.mDelta(ix(end))];
    %
    % Biomass spectrum:
    %
    %nexttile
    BB = B(:,ix-sim.p.idxB+1);
    BB = [BB(1,:); BB]; % Add dummy layer on top
    BB( BB<0.01 ) = 0.01;
    %panelField(m, -z, (B(:,ix-sim.p.idxB+1))');
    contourf( sim.p.m(sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup)), -z, BB, ...
        10.^linspace(-2,2,20),'linestyle','none')
    
    set(gca,'xscale','log','colorscale','log')    
    set(gca,'xtick',10.^(-9:2), 'XTickLabel',[])
    
    %caxis([0.01 100])
    %caxis([-3 2])
    
    title(sim.p.nameGroup(iGroup))
    if (iGroup==1)
        ylabel('Depth (m)')
    else
        set(gca,'yticklabel',[])
    end
    ylim(ylimit)
    
    if (iGroup == sim.p.nGroups)
        cbar = colorbar;
        cbar.Label.String  = 'biomass (\mug C l^{-1})';
        %set(cbar,'limits',[0.01, 100], ...
        %'ticks',[0.01,0.1,1,10,100],'ticklabels',{'0.01','0.1','1','10','100'})
    end
    %
    % Trophic strategy:
    %
    nexttile
    
    panelField(m,-z,colStrategy(ix-sim.p.idxB+1,:,:));
    set(gca,'xscale','log')
    set(gca,'xtick',10.^(-9:2:2))
    xlabel('Cell mass (\mugC)')
    if (iGroup==1)
        ylabel('Depth (m)')
    else
        set(gca,'yticklabel',[])
    end
    ylim(ylimit)
    
end
