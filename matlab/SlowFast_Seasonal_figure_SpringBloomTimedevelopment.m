% This figure plots one seasonal simulation with a specific mortHTL
%% parameters
%choose which mortHTL to use:
whichyear=[11.12:0.02:11.36];
% t=tiledlayout(length(whichyear),1,'TileIndexing','columnmajor');
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
mortHTL=0.5;

%% load simulation
load(fullfile('ModelResults','seasonal',['workspace_mortHTL_0_',num2str(10*mortHTL),'.mat']));
%% Plot figure
cm = parula(length(sim.p.nameGroup)-1);
marklines=(0.4*1d6*1d-12).*(4/3).*pi.*[2 20].^(3);
% figure;
% subplot(3,1,1)
% nexttile


% %% figure from one time
% 
% figure('Color','w');
% t=tiledlayout(length(whichyear),1,'TileIndexing','columnmajor');
% t.TileSpacing = 'tight';
% t.Padding = 'compact';
% 
% for ii=1:length(whichyear)
%     [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
%     themax(ii)=max(B_phag(plotthistime,:));
% end
% maxbio=max(themax);
% for ii=8 %1:length(whichyear)
%     nexttile
%     [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
% 
% 
%     % subplot(2,1,1)
%     hold on
%     mLower=p.mLower(p.ixStart(2):p.ixEnd(2));
%     mUpper=p.mUpper(p.ixStart(2):p.ixEnd(2));
%     remade=nan(length(sim.p.nameGroup)-1,n);
%     for i=2:length(sim.p.nameGroup)
%         remade(i-1,:)=B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2);
%     end
%     [maxRemade,maxRemadeIDX]=max(remade);
%     for i=1:length(mUpper)
%         mrange=[mLower(i) mUpper(i)];
%         fill([mrange fliplr(mrange)],[maxRemade(i)  maxRemade])
%     end
% 
% 
% 
% 
% 
%     for i=length(sim.p.nameGroup):-1:2
%         % plot(p.m(p.ixStart(i):p.ixEnd(i)),B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2),'Color',cm(i-1,:),'LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
%         area(p.m(p.ixStart(i):p.ixEnd(i)),B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2),'LineStyle','none','FaceColor',cm(i-1,:),'FaceAlpha',1,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
%         set(gca,'XScale','log');
%     end
% 
% 
%     xlim([min(p.mLower(p.ixStart(i):p.ixEnd(i))) max(p.mUpper(p.ixStart(i):p.ixEnd(i)))])
%     % theselims=ylim;
%     theselims=[0 maxbio];
%     % ylim([0 max(theselims)])
%     ylim([0 122])
%     for i=1:length(theselims)
%         plot([marklines(i) marklines(i)],[0 max(theselims)],'--','Color','r')
%     end
%     % legend('show','Location','eastoutside')
%     if mortHTL==0.1
%     ylabel(['t=',num2str(whichyear(ii))])
%     end
%     ax=gca;
%     if ii==length(whichyear)
%         xlabel(['mortHTL=',num2str(mortHTL)])
%     else
%         ax.XTickLabel=' ';
%     end
%     ax.Layer = 'top';
%     % ax.XGrid='on';
%     % ax.LineWidth=1.5;
%     % ax.FontSize=12;
% 
% end
% 
% 
% sgtitle(['B_{phag} spring bloom for at time t=',num2str(min(whichyear)),' to ',num2str(max(whichyear))])
% 

%%
Bphagidx_timebloom=nan(length(whichyear),n);
Bphag_timebloom=Bphagidx_timebloom;
for ii=1:length(whichyear)
    [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
    remade=nan(length(sim.p.nameGroup)-1,n);

    for i=2:length(sim.p.nameGroup)
        remade(i-1,:)=B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2);
    end
    [remadeMax_val,naxindex]=max(remade);
    Bphagidx_timebloom(ii,:)=naxindex;
    Bphag_timebloom(ii,:)=remadeMax_val;
end
Bphagidx_timebloom(Bphag_timebloom<0.1)=NaN;
figure;
pcolor(Bphagidx_timebloom')
%%
figure;
hold on
for ii=1:length(whichyear)
    [~,plotthistime]=min(abs(sim.t(ix)-(whichyear(ii).*365)));
    for i=2:length(sim.p.nameGroup)
        % B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2);
        Vv(i-1,:)=B_phag(plotthistime,p.ixStart(i)-2:p.ixEnd(i)-2);
        Vv(Vv<0.1)=NaN;
    end
    for tt=1:20
        [zz_sort,zz_idx]=sort(Vv(:,tt),'descend');

        for ll=1:5
            scatter(ii,p.m(p.ixStart(2)+tt-1),zz_sort(ll).*5,MarkerFaceColor=cm(zz_idx(ll),:),MarkerEdgeColor="none")
        end
    end

end
set(gca,'YScale','log');






