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
panelGlobal(sim.x,sim.y,log10(field),[0 2],sTitle='Surface DOC (log_{10} {\mu}g_C/l)',sProjection=sProjection);
set(gca,'XTickLabel','')

% Nitrogen
nexttile
field = getfield(sim.N);
c = panelGlobal(sim.x,sim.y,log10(field),[0 2],sTitle='Surface N (log_{10} {\mu}g_C/l)',sProjection=sProjection);
c.Label.String  = ' [\mug N l^{-1}]';
set(gca,'XTickLabel','')

% Silicate
if bSilicate
    nexttile
    field = getfield(sim.Si);
    c = panelGlobal(sim.x,sim.y,log10(field),[0 2],sTitle='Surface Si (log_{10} {\mu}g_{Si}/l)',sProjection=sProjection);
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
        [-2 1], sTitle=strcat(sim.p.nameGroup(i)), sProjection=sProjection);
    clim([-2 1])
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
        sTitle='Average net primary production (log_{10} g_C/m^2/yr)', ...
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



end


