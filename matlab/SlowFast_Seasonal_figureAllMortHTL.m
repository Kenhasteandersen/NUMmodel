% This figure plots one seasonal simulation with a specific mortHTL
%% parameters
%choose which mortHTL to use:
whichyear=[11.12:0.02:11.36];
figure('Color','w')
tiledlayout(length(whichyear),length(0.1:0.1:0.8),'TileIndexing','columnmajor');
for mortHTL=0.1:0.1:0.8

%% load simulation
load(fullfile('ModelResults','seasonal',['workspace_mortHTL_0_',num2str(10*mortHTL),'.mat']));
%% Plot figure
cm = parula(length(sim.p.nameGroup)-1);
% figure;
% subplot(3,1,1)
% nexttile


%% figure from one time

% figure('Color','w');
% tiledlayout(length(whichyear),1)

for ii=1:length(whichyear)
    [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
    themax(ii)=max(B_phag(plotthistime,:));
end
maxbio=max(themax);
for ii=1:length(whichyear)
    nexttile
    [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
    marklines=(0.4*1d6*1d-12).*(4/3).*pi.*[2 20].^(3);

    % subplot(2,1,1)
    hold on
    for i=2:length(sim.p.nameGroup)
        plot(p.m(p.ixStart(i):p.ixEnd(i)),B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2),'Color',cm(i-1,:),'LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
        set(gca,'XScale','log');
    end
    xlim([min(p.mLower(p.ixStart(i):p.ixEnd(i))) max(p.mUpper(p.ixStart(i):p.ixEnd(i)))])
    % theselims=ylim;
    theselims=[0 maxbio];
    % ylim([0 max(theselims)])
    ylim([0 122])
    for i=1:length(theselims)
        plot([marklines(i) marklines(i)],[0 max(theselims)],'--','Color','r')
    end
    % legend('show','Location','eastoutside')
    if mortHTL==0.1
    ylabel(['t=',num2str(whichyear(ii))])
    end
    if ii==length(whichyear)
        xlabel(['mortHTL=',num2str(mortHTL)])
    end
end
end

sgtitle(['B_{phag} spring bloom for at time t=',num2str(min(whichyear)),' to ',num2str(max(whichyear))])
%% related rates
figure('Color','w');
% tiledlayout(length(whichyear),1)
plotrates={'jDOC' 'jLreal' 'jFreal' 'mortpred' 'mort2'};
tiledlayout(length(whichyear),length(plotrates),'TileIndexing','columnmajor')
for jj=1:length(plotrates)
    thedata=plotrates{jj};

for ii=1:length(whichyear)

    [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
    mysim=sim.rates(ix(plotthistime));
    whichdata=mysim.(thedata);

    themax(ii)=max(whichdata);
end
maxbio=max(themax);
for ii=1:length(whichyear)
    nexttile
    [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
    marklines=(0.4*1d6*1d-12).*(4/3).*pi.*[2 20].^(3);
    mysim=sim.rates(ix(plotthistime));
    whichdata=mysim.(thedata);

    % subplot(2,1,1)
    hold on
    for i=2:length(sim.p.nameGroup)
        plot(p.m(p.ixStart(i):p.ixEnd(i)),whichdata(p.ixStart(i)-2:p.ixEnd(i)-2),'Color',cm(i-1,:),'LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
    end
    set(gca,'XScale','log');
    if strcmp(thedata,'mortpred')
        set(gca,'YScale','log');
    end

    xlim([min(p.mLower(p.ixStart(i):p.ixEnd(i))) max(p.mUpper(p.ixStart(i):p.ixEnd(i)))])
    % theselims=ylim;
    theselims=[0 maxbio];
    ylim([0 max(theselims)])
    for i=1:length(theselims)
        plot([marklines(i) marklines(i)],[0 max(theselims)],'--','Color','r')
    end
    % legend('show','Location','eastoutside')
    if jj==1
    ylabel(['t=',num2str(whichyear(ii))])
    end
    if ii==length(whichyear)
        xlabel(thedata)
    end
end
end
sgtitle(['bloom for mortHTL = ',num2str(mortHTL),' at time t=',num2str(min(whichyear)),' to ',num2str(max(whichyear))])








