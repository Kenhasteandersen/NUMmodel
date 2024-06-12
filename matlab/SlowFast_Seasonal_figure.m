% This figure plots one seasonal simulation with a specific mortHTL
%% parameters
%choose which mortHTL to use:
mortHTL=0.1;

%% load simulation
load(fullfile('ModelResults','seasonal',['workspace_mortHTL_0_',num2str(10*mortHTL),'.mat']));
%% Plot figure
cm = parula(length(sim.p.nameGroup)-1);
figure;
subplot(3,1,1)
area(sim.t(ix)./365,sum(sim.B(ix,:),2),'LineStyle','none','FaceColor',[0.901960784313726 0.901960784313726 0.901960784313726],'DisplayName','total biomass')
hold on
for i=2:length(sim.p.nameGroup)

    plot(sim.t(ix)./365,sum(sim.B(ix,p.ixStart(i)-2:p.ixEnd(i)-2),2),'-','LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))],'Color',cm(i-1,:))
end
set(gca,'XGrid','on','GridAlpha',1,'GridLineWidth',1,'XColor','r','XMinorTick','on','XMinorGrid','on','MinorGridLineStyle','--','MinorGridAlpha',1,'MinorGridLineWidth',0.5)
ax=gca; ax.XAxis.MinorTickValues=min(ax.XTick)-1:0.25:max(ax.XTick);
legend('show','Location','eastoutside')
subtitle('total biomass')

subplot(3,1,2)

area(sim.t(ix)./365,sum(B_phag(:,:),2),'LineStyle','none','FaceColor',[0.901960784313726 0.901960784313726 0.901960784313726],'DisplayName','total phagotrophs')
hold on
for i=2:length(sim.p.nameGroup)
    plot(sim.t(ix)./365,sum(B_phag(:,p.ixStart(i)-2:p.ixEnd(i)-2),2),'-','LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))],'Color',cm(i-1,:))
end
set(gca,'XGrid','on','GridAlpha',1,'GridLineWidth',1,'XColor','r','XMinorTick','on','XMinorGrid','on','MinorGridLineStyle','--','MinorGridAlpha',1,'MinorGridLineWidth',0.5)
ax=gca; ax.XAxis.MinorTickValues=min(ax.XTick)-1:0.25:max(ax.XTick);
legend('show','Location','eastoutside')
subtitle('phagotrophs')

subplot(3,1,3)

area(sim.t(ix)./365,sum(B_photo(:,:),2),'LineStyle','none','FaceColor',[0.901960784313726 0.901960784313726 0.901960784313726],'DisplayName','total phototrophs')
hold on
subtitle('phototrophs')
for i=2:length(sim.p.nameGroup)
    plot(sim.t(ix)./365,sum(B_photo(:,p.ixStart(i)-2:p.ixEnd(i)-2),2),'-','LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))],'Color',cm(i-1,:))
    % legendvector=[legendvector, sim.p.nameGroup{i}];
end
set(gca,'XGrid','on','GridAlpha',1,'GridLineWidth',1,'XColor','r','XMinorTick','on','XMinorGrid','on','MinorGridLineStyle','--','MinorGridAlpha',1,'MinorGridLineWidth',0.5)
ax=gca; ax.XAxis.MinorTickValues=min(ax.XTick)-1:0.25:max(ax.XTick);
legend('show','Location','eastoutside')
subtitle('phototrophs')

sgtitle(['Results for mortHTL = ',num2str(mortHTL)])

%% figure from one time
% whichyear=[11.2:0.02:11.4];
% whichyear=[11:0.1:12];
whichyear=[11.0:0.04:11.36];

figure('Color','w');
tiledlayout(length(whichyear),1)

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
    ylim([0 max(theselims)])
    for i=1:length(theselims)
        plot([marklines(i) marklines(i)],[0 max(theselims)],'--','Color','r')
    end
    % legend('show','Location','eastoutside')
    ylabel(['t=',num2str(whichyear(ii))])
end
sgtitle(['bloom for mortHTL = ',num2str(mortHTL),' at time t=',num2str(min(whichyear)),' to ',num2str(max(whichyear))])
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








