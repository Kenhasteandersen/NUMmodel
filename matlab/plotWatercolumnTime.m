%
% Plot a water column as a function of time at latitude `lat`
% and longitude `lon`. If lat and lon are empty it is assumed that the
% simulation is from "simulationWatercolumn".
%
% In:
%  sim - simulation structure either from global or from watercolumn
%  lat, lon - latitude and longitude (only to extract from global run)
% Optional named arguments:
%  bNewplot - boolean that decides whether to setup the plot
%  depthMax - mx depth for ylimit.
%
function plotWatercolumnTime(sim, lat, lon, options)

arguments
    sim struct;
    lat double = [];
    lon double = [];
    options.bNewPlot logical = true
    options.depthMax {mustBePositive} = [];
    options.nLevels = 20;
end

switch sim.p.nameModel
    case 'global'
        % Extract the water column from a global run:
        idx = calcGlobalWatercolumn(lat,lon,sim);
        N = squeeze(double(sim.N(idx.x, idx.y, idx.z,:)));
        DOC = squeeze(double(sim.DOC(idx.x, idx.y, idx.z,:)));
        if isfield(sim,'Si')
            Si = squeeze(double(sim.Si(idx.x, idx.y, idx.z,:)));
        end
        for i = 1:sim.p.nGroups
            B(i,:,:) = squeeze((sum(double(sim.B(idx.x, idx.y, idx.z,...
                (sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1,:)),4)));
        end
        z = sim.z(idx.z);
    case 'watercolumn'
        N = sim.N;
        DOC = sim.DOC;
        if isfield(sim,'Si')
            Si = sim.Si;
        end
        for i = 1:sim.p.nGroups
            B(i,:,:) = squeeze(sum(sim.B(:,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1,:),2));
        end
        z = sim.z + 0.5*sim.dznom;
    otherwise
        fprintf('Model type %s not supported.\n', sim.p.nameModel) 
end
N(N<=0) = 1e-8;
DOC(DOC<=0) = 1e-8;

t = sim.t;
%
% Make a layer at z = 0 with the same value as in the first grid point:
%
z = [0; z];
N = [N(1,:); N];
if isfield(sim,'Si')
    Si = [Si(1,:); Si];
end
DOC = [DOC(1,:); DOC];
B(:,2:length(z),:) = B;
B(:,1,:) = B(:,2,:);



if options.bNewPlot
    clf
    tiledlayout(2+isfield(sim,'Si')+sim.p.nGroups,1,'tilespacing','compact','padding','compact')
end

if isempty(options.depthMax)
    options.depthMax = max(z);
end
ylimit = [-options.depthMax, 0];

nexttile
%z = [sim.z-0.5*sim.dznom; sim.z(end)+0.5*sim.dznom(end)];
%panelField([t t(end)], -z, N');
%surface(t,-z, N)
contourf(t,-z,log10(N),linspace(-2,3,options.nLevels),'LineStyle','none')
title(['Nitrogen, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%set(gca,'ColorScale','log')
%shading interp
axis tight
colorbar('ticks',-2:3)
%caxis([-1 2])
ylim(ylimit)

if isfield(sim,'Si')
    nexttile
    %surface(t,-z, Si)
    contourf(t,-z,log10(Si),linspace(-2,3,options.nLevels),'LineStyle','none')
    title(['Silicate, lat ', num2str(lat),', lon ', num2str(lon)])
    ylabel('Depth (m)')
    %set(gca,'ColorScale','log')
    %shading interp
    axis tight
    colorbar('ticks',-2:3)
    %caxis([0.1 1000])
    ylim(ylimit)
end

nexttile
contourf(t,-z,log10(DOC),options.nLevels,'LineStyle','none')
%surface(t,-z, DOC)
title(['DOC, lat ', num2str(lat),', lon ', num2str(lon)])
ylabel('Depth (m)')
%xlabel('Concentration (\mugC l^{-1})')
%set(gca,'colorscale','log')
%set(gca,'ColorScale','log')
shading interp
axis tight
colorbar
%caxis([0.1,2])
ylim(ylimit)

for i = 1:sim.p.nGroups
    nexttile
    %surface(t,-z, squeeze(B(i,:,:)))
    B(B < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
    contourf(t,-z,log10(squeeze(B(i,:,:))),[linspace(-2,3,options.nLevels)],'LineStyle','none')
    title(sim.p.nameGroup(i))
    ylabel('Depth (m)')
    xlabel('Time (days)')
    %set(gca,'ColorScale','log')
    %shading interp
    axis tight
    colorbar
    %caxis([0.1 100])
    colorbar('ticks',-2:3,'limits',[-2 3])
    ylim(ylimit)
end

end