%
% Plot maps of global run
%
% In:
%  sim: the simulation to plot
%  iTime: The time step to plot. If iTime=0 (default) then it
%         plots the average.
%  sProjection: default 'fast'.
%
function plotGlobal(sim, iTime, options)

arguments
    sim struct;
    iTime double = 0;
    options.sProjection string = 'fast';
end
sProjection = options.sProjection;

if iTime == 0
    bAverageTime = true;
else
    bAverageTime = false;
end
%
% Do the plots:
%
clf
set(gcf,'color','w');

bSilicate = isfield(sim.p,'idxSi');

tiledlayout(2+bSilicate+sim.p.nGroups,1)
% DOC
%text(0, 1, labels(i),'Units','normalized')
nexttile
field = getfield(sim.DOC);
panelGlobal(sim.x,sim.y,log10(field),[-2 0],sTitle='Surface DOC',sProjection=sProjection);
set(gca,'XTickLabel','')

% Nitrogen
nexttile
field = getfield(sim.N);
c = panelGlobal(sim.x,sim.y,log10(field),[-2 1],sTitle='Surface N',sProjection=sProjection);
c.Label.String  = ' [\mug N l^{-1}]';
set(gca,'XTickLabel','')

% Silicate
if bSilicate
    nexttile
    field = getfield(sim.Si);
    c = panelGlobal(sim.x,sim.y,log10(field),[-2 1],sTitle='Surface Si',sProjection=sProjection);
    c.Label.String  = ' [\mug Si l^{-1}]';
end
set(gca,'XTickLabel','')

% Groups
for i = 1:sim.p.nGroups
    nexttile
    field = squeeze(sum(sim.B(:,:,:,:,(sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1),5));
    field = calcIntegrateGlobal(sim,field,bAverageTime);
    if ~bAverageTime
        field = squeeze(field(iTime,:,:));
    end
    cbar = panelGlobal(sim.x,sim.y,log10(field),...
        [0 1], sTitle=strcat('Surface log10(',sim.p.nameGroup(i),')'), sProjection=sProjection);
    clim([0 1])
    cbar.Label.String  = 'g_C/m^2';
    set(gca,'XTickLabel','')

    % Multicellular plankton
    %subplot(nPanels,1,nPanels)
    %panelGlobal(sim.x,sim.y,log10(sum(sim.B(:,:,1,findIxMulticellular(sim.p),iTime),4)),'Multicellular plankton (log10)',sProjection);
    %caxis([1 3])
end

if isfield(sim,'CnetPerArea')
    subplot(4,1,4)
    panelGlobal(sim.x, sim.y, log10(sim.CnetPerArea(:,:,1)), [0 3],...
        sTitle='Average net primary production (log10 gC/m2/yr)', ...
        sProjection=sProjection);
    clim([8 11])
end
set(gca,'xticklabel','auto')

    function field = getfield(fld)

        if bAverageTime
            field = squeeze(mean(fld(:,:,:,1),1));
        else
            field = squeeze(fld(iTime,:,:,1));
        end
        field(field<0) = 1e-20;
    end

%
% Plot transect
%
Lat = [45.01 41.93 38.63 35.3 31.55 28.61 25.05 21.45 17.74 13.81 10.1 6.26 -1.15 -4.985 -8.765 -12.085 -15.035 -21.6 -24.905 -27.2 -30.01 -32.375 -34.58 -37.01 -39.21 -41.37 -43.79 -46.02];
Lon = [-13.58 -16.02 -18.53 -20.92 -22.42 -24.36 -26.05 -27.75 -28.93 -28.29 -27.29 -26.38 -24.99 -24.97 -24.95 -24.93 -25.03 -25.05 -25.903 -28.25 -31.15 -33.67 -36.09 -38.83 -41.47 -44.01 -46.95 -49.87];
figure
plotGlobalTransect(sim,Lat,Lon)
sgtitle('Approximate AMT track - average over 1 year')

    end


