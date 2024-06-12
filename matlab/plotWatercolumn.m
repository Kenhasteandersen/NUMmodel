% θ
% Plot a water column from either a global or a water column simulation.
%
% Warning: the "setup" function needs to be called before this plot is
% made. If not matlab may crash or the results be incorrect.
%
% All biomasses are normalized to "Sheldon" spectra by dividing biomass
% with ln(Delta) (the ratio between the upper and lower mass in each mass
% bin); see Andersen and Visser (2023) Box V.
%
% In:
%  sim
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
    B = squeeze(double(sim.B(iTime, idx.x, idx.y, idx.z, :)));
    for i = 1:length(idx.z)
        if isfield(sim,'Si')
            u(i,:) = [sim.N(iTime, idx.x, idx.y, idx.z(i)), ...
                sim.DOC(iTime, idx.x, idx.y, idx.z(i)), ...
                sim.Si(iTime, idx.x, idx.y, idx.z(i)), ...
                B(i,:)];
        else
            u(i,:) = [sim.N(iTime, idx.x, idx.y, idx.z(i)), ...
                sim.DOC(iTime, idx.x, idx.y, idx.z(i)), B(i,:)];
        end
        L(i) = sim.L(iTime, idx.x, idx.y, idx.z(i));
        T(i) = sim.T(iTime, idx.x, idx.y, idx.z(i));
    end

else
    % Extract data from water column simulation:
    z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
    B = squeeze(double(sim.B(iTime, :, :)));

    % ixAve = find( sim.t > sim.t(end)/2 );
    % B = mean(squeeze(double(sim.B(:, :, ixAve))),2);

    for i = 1:length(sim.z)
        if isfield(sim,'Si')
            u(i,:) = [sim.N(iTime, i), sim.DOC(iTime,i), ...
                sim.Si(iTime,i), B(i,:)];
        else
            u(i,:) = [sim.N(iTime,i), sim.DOC(iTime,i), B(i,:)];
        end
        L(i) = sim.L(iTime,i);
        T(i) = sim.T(iTime,i);
    end
end

B(B<0) = 0;
%
% Calclulate trophic strategy and feeding level:
%
colStrategy = zeros(length(u)-sim.p.idxB+1, length(z)-1, 3);
colFeeding = colStrategy;
for i = 1:length(z)-1
    rates = getRates(sim.p, u(i,:), L(i), T(i));
    for j=1:size(u,2)-sim.p.idxB+1
        % Trophic strategy:
        colStrategy(j,i,:) = ...
            [min(1, max(0, 6*(rates.jFreal(j)))), ...
            min(1, max(0, 3*(rates.jLreal(j)))), ...
            min(1, max(0, 3*(rates.jDOC(j))))];
        % Feeding level:
        f(j,i) = rates.jFreal(j)./rates.jFmaxx(j);
        if isnan(f(j,i))
            colFeeding(j,i,:) = [0 0 0];
        else
            fc = rates.jR/rates.jMax; % Critical feeding level
            if f(j,i) < fc
                colFeeding(j,i,:) = [0, 1, f(j,i)/fc]; % Below critical feeding level
            else
                colFeeding(j,i,:) = [f(j,i), 0, 0]; % Above critical feeding level
            end
        end
    end
end
%
% Setup tiles:
%
if options.bNewplot
    clf
    tiles = tiledlayout(2,sim.p.nGroups,'tilespacing','compact','padding','compact');
    tiles.TileIndexing = 'columnmajor';
end

if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];
%
% Run over all groups.
%

Zmax = -inf;
Zmin = inf;

for iGroup = 1:sim.p.nGroups
    ix = sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup);
    m = [sim.p.mLower(ix), sim.p.mLower(ix(end))+sim.p.mDelta(ix(end))];
    %
    % Biomass spectrum:
    %
    if (~options.bNewplot) && (iGroup==1)
        h(iGroup) = gca;
    else
        h(iGroup) = nexttile;
    end
    %     BB = B(:,ix-sim.p.idxB+1);
    %     BB = [BB(1,:); BB]; % Add dummy layer on top
    %     BB( BB<0.01 ) = 0.01;

    % ix = sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup);
    % m = sim.p.m(ix);
    Delta = sim.p.mUpper(ix)./sim.p.mLower(ix);
    ixB = ix-sim.p.idxB+1;

    BB = B(:, ixB)./log(Delta); % Normalize to Sheldon spectra
    BB = [BB(1,:); BB]; % Add dummy layer on top
    BB( BB<0.01 ) = 0.01;

    % Treat groups with only one size class correctly:
    if size(BB,2)==1
        mm = [sim.p.mLower(sim.p.ixStart(iGroup)), sim.p.mUpper(sim.p.ixStart(iGroup))];
        BB(:,2) = BB;
    else
        mm = sim.p.m(sim.p.ixStart(iGroup):sim.p.ixEnd(iGroup));
    end
    contourf( mm, -z, BB, 10.^linspace(-2,2,20),'linestyle','none')

    set(gca,'xscale','log','colorscale','log')
    set(gca,'xtick',10.^(-9:2), 'XTickLabel',[])

    caxis([0.01 100])

    title(sim.p.nameGroup(iGroup))
    if (iGroup==1)
        ylabel('Depth (m)')
    else
        set(gca,'yticklabel',[])
    end
    ylim(ylimit)

    if max(max(BB(:,:))) > Zmax

        Zmax = max(max(BB(:,:)));
    end

    if min(min(BB(:,:))) < Zmin

        Zmin = min(min(BB(:,:)));
    end

    if (iGroup == sim.p.nGroups)
        cbar = colorbar;
        cbar.Label.String  = 'Sheldon biomass (\mug C l^{-1})';
        %         set(cbar,'limits',[0.01, 100], ...
        %         'ticks',[0.01,0.1,1,10,100],'ticklabels',{'0.01','0.1','1','10','100'})
        if Zmin==Zmax
            Zmax=Zmin+0.0001;
        end
        set(h, 'Colormap', jet, 'CLim', [Zmin Zmax])
    end
    %
    % Trophic strategy or feeding level:
    %
    nexttile

    if (sim.p.typeGroups(iGroup) < 10) % Unicellular
    panelField(m,-z,colStrategy(ix-sim.p.idxB+1,:,:));
    else
    panelField(m,-z,colFeeding(ix-sim.p.idxB+1,:,:));
    end
    set(gca,'xscale','log')
    if (sim.p.typeGroups(iGroup) < 10) | (sim.p.typeGroups(iGroup) >=100) % Unicellular
        set(gca,'xtick',10.^(-9:2:5))
        xlabel('Cell size (\mugC)')
    else
        set(gca,'xtick',10.^(-9:1:5))
        xlabel('Body size (\mugC)')
    end
    
    if (iGroup==1)
        ylabel('Depth (m)')
    else
        set(gca,'yticklabel',[])
    end
    ylim(ylimit)

end


if strcmp(sim.p.nameModel, 'watercolumn')
    sgtitle(['Sheldon biomass ({\mu}gC/l) at day: ', num2str(time), ', lat = ', num2str(sim.lat), char(176), ', lon = ', num2str(sim.lon), char(176)])
end
