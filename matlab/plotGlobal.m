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
panelGlobal(sim.x,sim.y,field,10.^linspace(0,2,10),...
    sTitle='Surface DOC',...
    sUnits='{\mu}g_C/l', sProjection=sProjection);
set(gca,'XTickLabel','','colorscale','log')


% Nitrogen
nexttile
field = getfield(sim.N);
panelGlobal(sim.x,sim.y,field, 10.^linspace(0,2,10), ...
    sTitle='Surface N', sUnits='{\mu}g_N/l', sProjection=sProjection);
set(gca,'XTickLabel','','colorscale','log')

% Silicate
if bSilicate
    nexttile
    field = getfield(sim.Si);
    panelGlobal(sim.x,sim.y,field,10.^linspace(0,2,10),...
        sTitle='Surface Si', sUnits='{\mu}g_{Si}/l', sProjection=sProjection);
    set(gca,'XTickLabel','','colorscale','log')
end

% Groups
for i = 1:sim.p.nGroups
    nexttile
    field = squeeze(sum(sim.B(:,:,:,:,(sim.p.ixStart(i):sim.p.ixEnd(i))-sim.p.idxB+1),5));
    field = calcIntegrateGlobal(sim,field,bAverageTime);
    if ~bAverageTime
        field = squeeze(field(iTime,:,:));
    end
    cbar = panelGlobal(sim.x,sim.y,field,...
        10.^linspace(-2, 1,10), sUnits='g_C/m^2',...
        sTitle=strcat(sim.p.nameGroup(i)), sProjection=sProjection);
    clim(10.^[-2 1])
    set(gca,'XTickLabel','','colorscale','log')

    % Multicellular plankton
    %subplot(nPanels,1,nPanels)
    %panelGlobal(sim.x,sim.y,log10(sum(sim.B(:,:,1,findIxMulticellular(sim.p),iTime),4)),'Multicellular plankton (log10)',sProjection);
    %caxis([1 3])
end

if isfield(sim,'CnetPerArea')
    subplot(4,1,4)
    panelGlobal(sim.x, sim.y, sim.CnetPerArea(:,:,1), 10.^linspace(0,3,10),...
        sTitle='Average net primary production', ...
        sProjection=sProjection, sUnits='g_C/m^2/yr');
    clim(10.^[8 11])
end
set(gca,'xticklabel','auto','colorscale','log')

    function field = getfield(fld)

        if bAverageTime
            field = squeeze(mean(fld(:,:,:,1),1));
        else
            field = squeeze(fld(iTime,:,:,1));
        end
        field(field<0) = 1e-20;
    end



end


