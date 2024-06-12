% This figure plots seasonal changes in biomass (total, phago and photo) for a specific 
% mortHTL. Data is calculated in SlowFast_seasonal_mortHTL_Loop.m
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