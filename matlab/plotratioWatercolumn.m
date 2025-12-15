% 
% Plot a water column as a function of time at latitude `lat`
% and longitude `lon`. If lat and lon are empty it is assumed that the
% simulation is from "simulationWatercolumn".
%
% In:
%  sim - simulation structure either from global or from watercolumn
%  lat, lon - latitude and longitude (only to extract from global run)
%
% Optional named arguments:
%  bNewplot - boolean that decides whether to setup the plot
%  depthMax - max depth.
%  bOnlyLastYear - boolean deciding whether to only plot the last year
%  bMolarUnits - boolean deciding whether to show nutrients in molar units.
%
function plotratioWatercolumn(sim, lat, lon, options)

arguments
    sim struct;
    lat double = 60;
    lon double = -40;
    options.bNewPlot logical = true;
    options.depthMax {mustBePositive} = [];
    options.bOnlyLastYear = false;
    options.nLevels = 20;
    options.bMolarUnits = true;
end

t = sim.t;

switch sim.p.nameModel
    case 'global'
        idx = calcGlobalWatercolumn(lat,lon,sim);
        for i = 1:sim.p.nGroups
            B(i,:,:) = squeeze((sum(double(sim.B(:,idx.x, idx.y, idx.z,...
                (sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1)),5)))';
        end
        iGel  = 3:8;
        iNon  = 9:14;
        Bgel  = squeeze(sum(B(iGel,:,:), 1));     % sum over groups 1–7
        Bnon  = squeeze(sum(B(iNon,:,:), 1));     % sum over groups 8–14
        RatioGel = Bgel ./ (Bnon+Bgel);   


        z = sim.z(idx.z);
    case 'watercolumn'
        for i = 1:sim.p.nGroups
            B(i,:,:) = squeeze(sum(sim.B(:,:,(sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1),3))';
        end
        iGel  = 3:7;
        iNon  = 8:12;
        Bgel  = squeeze(sum(B(iGel,:,:), 1));     % sum over groups 1–7
        Bnon  = squeeze(sum(B(iNon,:,:), 1));     % sum over groups 8–14
        RatioGel = Bgel ./ (Bnon+Bgel);   

        z = sim.z + 0.5*sim.dznom;
        lat = sim.lat;
        lon = sim.lon;
    otherwise
        fprintf('Model type %s not supported.\n', sim.p.nameModel)
end

%
% Make a layer at z = 0 with the same value as in the first grid point:
%
z = [0; z];
RatioGel = [RatioGel(1,:); RatioGel];

nTiles = 2;

if options.bNewPlot
    clf
    tiledlayout(nTiles,1,'tilespacing','tight','padding','tight')
end

if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];
xlimit = [sim.t(1) sim.t(end)];
if options.bOnlyLastYear
    xlimit(1) = sim.t(end)-365;
end

nexttile
contourf(t, -z, RatioGel, linspace(0,1,options.nLevels), 'LineStyle','none')
title('RatioGel','FontWeight','normal')
axis tight
h = colorbar;
h.Label.String = 'ratio';
caxis([0 1])
ylim(ylimit)
xlim(xlimit)
set(gca, 'colorscale','linear')   % important : échelle linéaire
set(gca,'XTickLabel','')

RatioGel_daily = mean(RatioGel, 1);   % 1 x nDays

% Dates et mois
dates = datetime(2024,1,1) + days(t-1);  
mois = month(dates);

mask3 = year(dates) >= year(dates(1)) + 2;   % Commence à l'année 3
dates_3 = dates(mask3);
RatioGel_daily_3 = RatioGel_daily(mask3);
mois_3 = month(dates_3);

% Boxplot mensuel
nexttile
boxplot(RatioGel_daily_3, mois_3)
xticklabels(month(datetime(2024,1:12,1),'name'))
ylabel('RatioGel relatif')
xlabel('Mois')
title('Distribution mensuelle de RatioGel')


if strcmp(sim.p.nameModel, 'watercolumn')
    sgtitle(['Water column at lat = ', num2str(lat), char(176), ', lon = ', num2str(lon), char(176)])
end
annotation('textbox', [0.075, 0.5, 0.5, 0.04], 'String', 'Depth (m)', 'FontSize', 10,'rotation',90,...
    'edgecolor','none','VerticalAlignment','bottom');

end

