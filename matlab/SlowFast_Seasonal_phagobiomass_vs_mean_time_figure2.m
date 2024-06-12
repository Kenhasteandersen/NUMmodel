clear;
% This figure plots seasonal changes in biomass (total, phago and photo) for a specific
% mortHTL. Data is calculated in SlowFast_seasonal_mortHTL_Loop.m
% ** for article
%% parameters
%choose which mortHTL to use:
mortHTL=0.5;

% min limits for values
minval=0.01;

% color for the pico-nano-micro grid
gridcolor=[0.8 0.8 0.8];
% choose which timesteps to use:
whichyear=linspace(49,131,8);
marklines=(0.4*1d6*1d-12).*(4/3).*pi.*[2 20].^(3);

%% load simulation
load(fullfile('ModelResults','seasonal',['workspace_mortHTL_0_',num2str(10*mortHTL),'.mat']));
% set  small values to minval
sim.B(sim.B<0.01)=minval;

%% Plot the different alpha and palatability groups
figure('Color','w');
cm = parula(length(sim.p.nameGroup)-1);
[alphaF,palatability]=getPalatability;
%incrase the resolution of alppaF and palatability
alphaF_highres=exp(interp(log(alphaF),100));
palatability_highres=interp1(alphaF,palatability,alphaF_highres);
plot(alphaF_highres,palatability_highres,':','Color',[0.8 0.8 0.8],'LineWidth',1.5)
hold on
scatter(alphaF_these,palatability_these,60,cm,'filled')
set(gca,'XScale','log');
xlabel('Specific clerance rate (l day^{-1} \mugC^{-1})')
ylabel('Vulnerability')
box off
ax=gca;
ax.LineWidth=1.5;
ax.FontSize=12;

%% Divide biomass into osmo, photo and phago
% set  small values to minval
sim.B(sim.B<0.01)=minval;
% divide biomass into osmo; photo and phago:
B_osmo=nan(size(sim.B));
B_photo=B_osmo;
B_phag=B_osmo;
for i=1:length(sim.B)
    iBiomass=sim.B(i,:);
    iDOC=sim.rates(i).jDOC';
    iF=sim.rates(i).jFreal';
    iL=sim.rates(i).jLreal';
    B_osmo(i,:)=iBiomass.*iDOC./(iDOC+iL+iF);
    B_phag(i,:)=iBiomass.*iF./(iDOC+iL+iF);
    B_photo(i,:)=iBiomass.*iL./(iDOC+iL+iF);

end
%% move biomass to a linear timeline
i=1;
t=linspace(sim.t(ix(i)),sim.t(ix(end)),5000);

B_interped=nan(5000,size(sim.B,2));
B_phag_interped=B_interped;
B_photo_interped=B_interped;
B_osmo_interped=B_interped;
for i=1:size(sim.B,2) % for all timesteps
    B_interped(:,i)=interp1(sim.t(ix),sim.B(ix,i),t);
    B_phag_interped(:,i)=interp1(sim.t(ix),B_phag(ix,i),t);
    B_photo_interped(:,i)=interp1(sim.t(ix),B_photo(ix,i),t);
    B_osmo_interped(:,i)=interp1(sim.t(ix),B_osmo(ix,i),t);
end

%% Plot biomass across the year
t_mod=floor(mod(t-1,365)+1);
BiomassAcrossTheYear(B_interped,sim,p,alphaF_these,t_mod,'Total biomass across year')
% BiomassAcrossTheYear(B_osmo_interped,sim,p,alphaF_these,t_mod,'osmotrophic biomass across year')
% BiomassAcrossTheYear(B_photo_interped,sim,p,alphaF_these,t_mod,'phototrophic biomass across year')
% BiomassAcrossTheYear(B_phag_interped,sim,p,alphaF_these,t_mod,'phagotrophic biomass across year')




%%
for ii=2:length(sim.p.nameGroup)
    for itime=1:365
        this_time=find(t_mod==itime);
        theseB=B_interped(this_time,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        % % take the mean of all dates
        B_meanYaar(itime,1:20,ii)=mean(theseB,1,'omitnan');
        B_minYaar(itime,1:20,ii)=min(theseB);
        B_maxYaar(itime,1:20,ii)=max(theseB);



        theseB_phag=B_phag_interped(this_time,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        % % take the mean of all dates
        B_phag_meanYaar(itime,1:20,ii)=mean(theseB_phag,1,'omitnan');
        B_phag_minYaar(itime,1:20,ii)=min(theseB_phag);
        B_phag_maxYaar(itime,1:20,ii)=max(theseB_phag);

        theseB_photo=B_photo_interped(this_time,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        B_photo_meanYaar(itime,1:20,ii)=mean(theseB_photo,1,'omitnan');
        B_photo_minYaar(itime,1:20,ii)=min(theseB_photo);
        B_photo_maxYaar(itime,1:20,ii)=max(theseB_photo);

        theseB_osmo=B_osmo_interped(this_time,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        B_osmo_meanYaar(itime,1:20,ii)=mean(theseB_osmo,1,'omitnan');
        B_osmo_minYaar(itime,1:20,ii)=min(theseB_osmo);
        B_osmo_maxYaar(itime,1:20,ii)=max(theseB_osmo);
    end
end
%% Plot biomass for spring
BiomassAcrossSpring(B_interped,whichyear,sim,p,marklines,alphaF_these,t_mod,'Spring bloom for total biomass',B_phag_meanYaar,B_photo_meanYaar,B_osmo_meanYaar,B_meanYaar)
% BiomassAcrossSpring(B_photo_interped,whichyear,sim,p,marklines,alphaF_these,t_mod,'Spring bloom for phototroph biomass')
%%
figure('color','w')
t=tiledlayout(length(whichyear),1,'TileIndexing','columnmajor');
t.TileSpacing ="loose";  
t.Padding = 'compact';
for ii=1:length(whichyear)
    nexttile(ii)
    thisyear=whichyear(ii);
    % subplot(length(whichyear),1,ii);
    xlim([min(p.mLower(p.ixStart(2):p.ixEnd(2))) max(p.mUpper(p.ixStart(2):p.ixEnd(2)))])
    ylim([ii-0.5 ii+0.5]);
    grid on
    set(gca,'XScale','log');
    ax=gca;
    if ii<length(whichyear)
    ax.XTickLabel={};
    end
    set(gca,'ytick',[])
    ax.YColor='w';
end
%%
figure('color','w')
t=tiledlayout(length(whichyear),length(sim.p.nameGroup)-1,'TileIndexing','columnmajor');
t.TileSpacing ="loose";
t.Padding = 'compact';
tnr=0;
for i=2:length(sim.p.nameGroup)%2:length(sim.p.nameGroup)
    for ii=1:length(whichyear)
        tnr=tnr+1;
        nexttile(tnr)
        thisyear=whichyear(ii);
        % subplot(length(whichyear),1,ii);
        xlim([min(p.mLower(p.ixStart(2):p.ixEnd(2))) max(p.mUpper(p.ixStart(2):p.ixEnd(2)))])
        ylim([ii-0.5 ii+0.5]);
        grid on
        set(gca,'XScale','log');
        ax=gca;
        if ii<length(whichyear)
            ax.XTickLabel={};
        end
        set(gca,'ytick',[])
        ax.YColor='w';

        biophag=squeeze(B_phag_meanYaar(round(thisyear),:,i))';

        biophot=squeeze(B_photo_meanYaar(round(thisyear),:,i))';
        bioosmo=squeeze(B_osmo_meanYaar(round(thisyear),:,i))';
        xx=(1-((biophot+bioosmo)./(biophag+biophot+bioosmo)));
        bio=squeeze(B_meanYaar(round(thisyear),:,i))';
        [bio]=removeTooSmallPoints(bio,0.5);
        xx(isnan(bio))=NaN;
        these=find(~isnan(xx));
        if ~isempty(these)
            for iii=1:length(these)
                mlower=p.mLower(p.ixStart(2):p.ixEnd(2));
                mUpper=p.mUpper(p.ixStart(2):p.ixEnd(2));
                pp=patch([mlower(these(iii)) mlower(these(iii)) mUpper(these(iii)) mUpper(these(iii))],[ii-0.5 ii+0.5 ii+0.5 ii-0.5],cm(i-1,:),'EdgeColor','none','FaceAlpha',xx(these(iii)),'AlphaDataMapping','scaled');
                % plot([mlower(these(iii)) mlower(these(iii)) mUpper(these(iii)) mUpper(these(iii))],[ii-0.5 ii+0.5 ii+0.5 ii-0.5],'FaceColor',cm(i-1,:),'AlphaDataMapping','scaled','AlphaData',xx(these))
                text((mUpper(these(iii))+mlower(these(iii)))./2,ii,num2str(xx(these(iii)).*100),'Rotation',90)
            end
        end
    end
end
sgtitle('fraction phago')

%% i pct i graf
figure('color','w')
t=tiledlayout(length(whichyear),1,'TileIndexing','columnmajor');
t.TileSpacing ="loose";
t.Padding = 'compact';
for i=2:length(sim.p.nameGroup)%2:length(sim.p.nameGroup)
    for ii=1:length(whichyear)
        nexttile(ii)
        thisyear=whichyear(ii);
        % subplot(length(whichyear),1,ii);
        xlim([min(p.mLower(p.ixStart(2):p.ixEnd(2))) max(p.mUpper(p.ixStart(2):p.ixEnd(2)))])
        ylim([ii-0.5 ii+0.5]);
        grid on
        set(gca,'XScale','log');
        % ax=gca;
        % if ii<length(whichyear)
        %     ax.XTickLabel={};
        % end
        % set(gca,'ytick',[])
        % ax.YColor='w';
        ax=gca; ax.XGrid='on';ax.XMinorGrid='off';
        ax.XTick=p.mUpper(p.ixStart(2):p.ixEnd(2));

        biophag=squeeze(B_phag_meanYaar(round(thisyear),:,i))';

        biophot=squeeze(B_photo_meanYaar(round(thisyear),:,i))';
        bioosmo=squeeze(B_osmo_meanYaar(round(thisyear),:,i))';
        xx=(1-((biophot+bioosmo)./(biophag+biophot+bioosmo))).*100;
        % bio=squeeze(B_meanYaar(round(thisyear),:,i))';
        % [bio]=removeTooSmallPoints(bio,0.5);
        % xx(isnan(bio))=NaN;
        these=find(~isnan(xx));
        m=p.m(p.ixStart(2):p.ixEnd(2));
        hold on;
        plot(m,xx,'*','Color',cm(i-1,:))
        ylim([0 100])
        % if ~isempty(these)
        %     for iii=1:length(these)
        %         plot(xx(these(iii)))
        %         mlower=p.mLower(p.ixStart(2):p.ixEnd(2));
        %         mUpper=p.mUpper(p.ixStart(2):p.ixEnd(2));
        %         pp=patch([mlower(these(iii)) mlower(these(iii)) mUpper(these(iii)) mUpper(these(iii))],[ii-0.5 ii+0.5 ii+0.5 ii-0.5],cm(i-1,:),'EdgeColor','none','FaceAlpha',xx(these(iii)),'AlphaDataMapping','scaled');
        %         % plot([mlower(these(iii)) mlower(these(iii)) mUpper(these(iii)) mUpper(these(iii))],[ii-0.5 ii+0.5 ii+0.5 ii-0.5],'FaceColor',cm(i-1,:),'AlphaDataMapping','scaled','AlphaData',xx(these))
        %         text((mUpper(these(iii))+mlower(these(iii)))./2,ii,num2str(xx(these(iii)).*100),'Rotation',90)
        %     end
        % end
    end
end
sgtitle('fraction phago')

%% Hvad er ratio mellem de tre?
for ii=1:length(whichyear)
    figure('Color','w','Name','ratio biomass')
    hold on
    thisyear=whichyear(ii);
    hold on
    for i=2:length(sim.p.nameGroup)
        biophag=squeeze(B_phag_meanYaar(round(thisyear),:,i))';
        nexttile
        biophot=squeeze(B_photo_meanYaar(round(thisyear),:,i))';
        bioosmo=squeeze(B_osmo_meanYaar(round(thisyear),:,i))';
        xx=((biophot+bioosmo)./(biophag+biophot+bioosmo));
        bio=squeeze(B_meanYaar(round(thisyear),:,i))';
        [bio]=removeTooSmallPoints(bio,0.5);
        xx(isnan(bio))=NaN;
        bar(xx,'FaceColor',cm(i-1,:));
        xlim([1 20])
        ax=gca;
        ax.XTick=1:20;
        ax.XTickLabel=p.m(p.ixStart(2):p.ixEnd(2));
        %     %
        ylim([0 1])
        title(['group nr.',num2str(i)])

    end
    sgtitle(['day ',num2str(round(whichyear(ii)))])
end
title('phagotrophs')


%% FUNCTIONS
function [Bio]=removeTooSmallPoints(Bio,minlim)

% if the first and the next index is too small, remove the first
if Bio(1:2)<minlim
    Bio(1)=NaN;
end
% if the index and the one before and after is too small (or nan), remove it
for jj=2:length(Bio)-1
    if Bio(jj)<minlim && Bio(jj+1)<minlim && (isnan(Bio(jj-1)) || Bio(jj-1)<minlim)
        Bio(jj)=NaN;
    end
end
if Bio(end)<minlim && isnan(Bio(end-1))
    Bio(end)=NaN;
end
end

function BiomassAcrossTheYear(B,sim,p,alphaF,t_mod,thistitle)
% Function that plots the biomass as a function of day in the year. It
% plots all the different alphaF groups on their own.
figure('Color','w')
cm = parula(length(sim.p.nameGroup)-1);
x2 = [1:365, fliplr(1:365)];
hold on
B_mean=nan(1,365);
B_min=B_mean;
B_max=B_mean;
for ii=2:length(sim.p.nameGroup)
    for itime=1:365
        this_time=t_mod==itime;
        theseB_phag=B(this_time,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        % sum all size classes and take the mean of all the dates:
        B_mean(itime)=mean(sum(theseB_phag,2,'omitnan'));
        % sum all size classes and take the minimum of all the dates:
        B_min(itime)=min(sum(theseB_phag,2,'omitnan'));
        % sum all size classes and take the maximum of all the dates:
        B_max(itime)=max(sum(theseB_phag,2,'omitnan'));
    end
    inBetween = [B_min, fliplr(B_max)];
    fill(x2, inBetween,'g','LineStyle','none','FaceColor',cm(ii-1,:),'FaceAlpha',0.2);
    plot(1:365,B_mean,'-','LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF(ii-1))],'Color',cm(ii-1,:))
end
xlabel('Days of year');
ylabel('Biomass (\mug C l^{-1})')
set(gca,'XGrid','on','GridAlpha',0.4,'GridLineWidth',0.5,'XMinorTick','on','XMinorGrid','on','MinorGridLineStyle',':','MinorGridAlpha',0.4,'MinorGridLineWidth',0.5)
theyear=cumsum([31 29 31 30 31 30 31 31 30 31 30 31]);
ax=gca;ax.XTick=[theyear(2) theyear(5) theyear(8) theyear(11)];
ax.XAxis.MinorTickValues=theyear;
xlim([0 365])
% general axis settings
ax.LineWidth=1.5;
ax.FontSize=12;
ax.Layer = 'top';
title(thistitle)
end

function BiomassAcrossSpring(B,whichyear,sim,p,marklines,alphaF_these,t_mod,titletext,B_phag_meanYaar,B_photo_meanYaar,B_osmo_meanYaar,B_meanYaar)
cm = parula(length(sim.p.nameGroup)-1);
B_mean=nan(365,20,length(sim.p.nameGroup));
B_min=B_mean;
B_max=B_mean;


for ii=2:length(sim.p.nameGroup)
    for i=1:365
        iix=t_mod==i;
        theseB_photo=B(iix,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        % get the mean, min and max of all "samples" from this day:
        B_mean(i,1:20,ii)=mean(theseB_photo,1,'omitnan');
        B_min(i,1:20,ii)=min(theseB_photo);
        B_max(i,1:20,ii)=max(theseB_photo);
    end
end

figure('Color','w','Name',titletext)
t=tiledlayout(length(whichyear),1,'TileIndexing','columnmajor');
t.TileSpacing = 'tight';
t.Padding = 'compact';

for ii=1:length(whichyear)
    tplace=0;
    thisyear=whichyear(ii);

    nexttile;hold on
    plot([marklines(1) marklines(1)],[0 100],'--','Color',[0.8 0.8 0.8])
    plot([marklines(2) marklines(2)],[0 100],'--','Color',[0.8 0.8 0.8])
   

    for i=2:length(sim.p.nameGroup)
        if i<5
            linetypes={'-','--'};
        else
            linetypes={':','-.-'};
        end
        % create grey area
        x_vector = [p.m(p.ixStart(i):p.ixEnd(i)), fliplr(p.m(p.ixStart(i):p.ixEnd(i)))];
        inBetween = [squeeze(B_min(round(thisyear),:,i)), fliplr(squeeze(B_max(round(thisyear),:,i)))];
        fill(x_vector', inBetween,'g','LineStyle','none','FaceColor',cm(i-1,:),'FaceAlpha',0.2);
        
        %remove the points with no biomass so it does not clutter the plot
        Bio=squeeze(B_mean(round(thisyear),:,i))';
        [Bio]=removeTooSmallPoints(Bio,0.5);
        TF = islocalmin(Bio);
        % ADD NUMBERS
        biophag=squeeze(B_phag_meanYaar(round(thisyear),:,i))';
        biophot=squeeze(B_photo_meanYaar(round(thisyear),:,i))';
        bioosmo=squeeze(B_osmo_meanYaar(round(thisyear),:,i))';
        xx=(1-((biophot+bioosmo)./(biophag+biophot+bioosmo)));
        bio=squeeze(B_meanYaar(round(thisyear),:,i))';
        [bio]=removeTooSmallPoints(bio,0.5);
        xx(isnan(bio))=NaN;
        these=find(~isnan(xx));

        m=p.m(p.ixStart(2):p.ixEnd(2));
        tplace=tplace+50;
        if sum(TF)==0

            plot(p.m(p.ixStart(i):p.ixEnd(i)),Bio,'Color',cm(i-1,:),'LineWidth',1.5,'LineStyle',linetypes{1},'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
            if ~isempty(these)
                for iii=1:length(these)
                    text(m(these(iii)),tplace,num2str(round(xx(these(iii)).*100)),'Color',cm(i-1,:))
                end
            end

        else
            difidx=find(TF==1,1,'last');
            plot(m(1:difidx),Bio(1:difidx),'Color',cm(i-1,:),'LineWidth',1.5,'LineStyle',linetypes{1},'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
            plot(m(difidx:end),Bio(difidx:end),'Color',cm(i-1,:),'LineWidth',1.5,'LineStyle',linetypes{2},'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
            % plot(m(TF),Bio(TF),'o','Color',cm(i-1,:),'MarkerSize',24)
            if ~isempty(these)
                for iii=1:length(these)
                    text(m(these(iii)),tplace,num2str(round(xx(these(iii)).*100)),'Color',cm(i-1,:))
                end
            end
        end

    end
    % wrap up axis and titel
    set(gca,'XScale','log');
    ax=gca;
    ax.Layer = 'top';
    ax.LineWidth=1;
    ax.FontSize=12;
    xlim([min(p.mLower(p.ixStart(2):p.ixEnd(2))) max(p.mUpper(p.ixStart(2):p.ixEnd(2)))])
    ylim([0 200])
    %for pct opmåling
    ax=gca; ax.XGrid='on';ax.XMinorGrid='off';
    ax.XTick=p.mUpper(p.ixStart(2):p.ixEnd(2));
    % ax.YAxis.TickLabels={'0';'';'100'};
    if ii<length(whichyear)
        ax=gca;
        ax.XAxis.TickLabels='';
    else
        xlabel('Cell mass (\mug C)')
    end
    ttt=title(['day ',num2str(round(whichyear(ii)))],'FontSize',12,'FontWeight','normal');
    ttt.HorizontalAlignment='right';
    ttt.Position=[p.m(p.ixEnd(2)) 10 0];
    ttt.FontAngle='italic';

end
end
