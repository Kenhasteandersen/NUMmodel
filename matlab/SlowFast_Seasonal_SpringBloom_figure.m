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


whichyear=linspace(49,131,8);

%% load simulation
load(fullfile('ModelResults','seasonal',['workspace_mortHTL_0_',num2str(10*mortHTL),'.mat']));
cm = parula(length(sim.p.nameGroup)-1);

%% set  small values to minval
sim.B(sim.B<minval)=minval;

%% divide into photo/osmo, mixo and phago
Bosmo_photo=nan(size(sim.B));
Bmixo=nan(size(sim.B));
Bphago=nan(size(sim.B));

for i=1:length(sim.B)
    iB=sim.B(i,:);
    iDOC=sim.rates(i).jDOC';
    iF=sim.rates(i).jFreal';
    iL=sim.rates(i).jLreal';
    

    iFratio=iF./(iDOC+iF+iL);
    iosmopotoratio=(iDOC+iL)./(iDOC+iF+iL);

    op_bry=0.01;
    mix_bry=0.6;

    % the places where phagotrophy is too small, register bio as photo+osmo 
    iBosmo_photo=iB;
    iBosmo_photo(iFratio>op_bry)=minval; %all other places use minval

    % the others are mixotrophic
    iBmixo=iB;
    iBmixo(iFratio<op_bry & iFratio>mix_bry)=minval; %all other places use minval

    % the places where phagotrophy is large, register bio as phago 
    iBphag=iB;
    iBphag(iFratio<mix_bry)=minval;%all other places use minval


    Bosmo_photo(i,:)=iBosmo_photo;
    Bmixo(i,:)=iBmixo;
    Bphago(i,:)=iBphag;


end

% sim.rates(1).jDOC
% B_phag(B_phag<minval)=minval;
% B_photo(B_photo<minval)=minval;

%% move biomass to a linear timeline
t=linspace(sim.t(ix(1)),sim.t(ix(end)),5000);
B_interped=nan(5000,size(sim.B,2));
Bosmo_photo_interped=B_interped;
Bphag_interped=B_interped;
Bmixo_interped=B_interped;
for i=1:size(sim.B,2) % for all timesteps
    B_interped(:,i)=interp1(sim.t(ix),sim.B(ix,i),t);
    Bphag_interped(:,i)=interp1(sim.t(ix),Bphago(ix,i),t);
    Bosmo_photo_interped(:,i)=interp1(sim.t(ix),Bosmo_photo(ix,i),t);
    Bmixo_interped(:,i)=interp1(sim.t(ix),Bmixo(ix,i),t);
end

t_mod=floor(mod(t-1,365)+1);

% %% figure of the alpha and palatability groups
% figure('Color','w');
% [alphaF,palatability]=getPalatability;
% %incrase the resolution of alppaF and palatability
% alphaF_highres=exp(interp(log(alphaF),100));
% palatability_highres=interp1(alphaF,palatability,alphaF_highres);
% plot(alphaF_highres,palatability_highres,':','Color',[0.8 0.8 0.8],'LineWidth',1.5)
% hold on
% scatter(alphaF_these,palatability_these,60,cm,'filled')
% set(gca,'XScale','log');
% xlabel('Specific clerance rate (l day^{-1} \mugC^{-1})')
% ylabel('Vulnerability')
% box off
% ax=gca;
% ax.LineWidth=1.5;
% ax.FontSize=12;
% % ax.Layer = 'top';



%% Get yearly mean across all years
B_interped=nan(5000,size(sim.B,2));
Bphag_mean=nan(365,20,length(sim.p.nameGroup));
Bphag_min=Bphag_mean;
Bphag_max=Bphag_mean;
Bosmo_photo_mean=Bphag_mean;
Bosmo_photo_min=Bphag_mean;
Bosmo_photo_max=Bphag_mean;
Bmixo_mean=Bphag_mean;
Bmixo_min=Bphag_mean;
Bmixo_max=Bphag_mean;

for ii=2:length(sim.p.nameGroup)
    for i=1:365
        iix=find(t_mod==i); %daily index
        iBphag=Bphag_interped(iix,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        iBosmo_photo=Bosmo_photo_interped(iix,p.ixStart(ii)-2:p.ixEnd(ii)-2);
        iBmixo=Bmixo_interped(iix,p.ixStart(ii)-2:p.ixEnd(ii)-2);

        Bphag_mean(i,1:20,ii)=mean(iBphag,1,'omitnan');
        Bphag_min(i,1:20,ii)=min(iBphag);
        Bphag_max(i,1:20,ii)=max(iBphag);

        Bosmo_photo_mean(i,1:20,ii)=mean(iBosmo_photo,1,'omitnan');
        Bosmo_photo_min(i,1:20,ii)=min(iBosmo_photo);
        Bosmo_photo_max(i,1:20,ii)=max(iBosmo_photo);

        Bmixo_mean(i,1:20,ii)=mean(iBmixo,1,'omitnan');
        Bmixo_min(i,1:20,ii)=min(iBmixo);
        Bmixo_max(i,1:20,ii)=max(iBmixo);
    end
end
%% Plot figure
figure('Color','w','Name','Spring bloom for phototroph biomass')
t=tiledlayout(length(whichyear),3,'TileIndexing','columnmajor');
t.TileSpacing = 'tight';
t.Padding = 'compact';
plotvertical(Bosmo_photo_min,Bosmo_photo_max,Bosmo_photo_mean,sim,p,cm,whichyear,alphaF_these)
plotvertical(Bmixo_min,Bmixo_max,Bmixo_mean,sim,p,cm,whichyear,alphaF_these)
plotvertical(Bphag_min,Bphag_max,Bphag_mean,sim,p,cm,whichyear,alphaF_these)
%%
function plotvertical(Bosmo_photo_min,Bosmo_photo_max,Bosmo_photo_mean,sim,p,cm,whichyear,alphaF_these)
marklines=(0.4*1d6*1d-12).*(4/3).*pi.*[2 20].^(3);
for ii=1:length(whichyear)
    nexttile
    hold on
    % plot nano and micro lines
    plot([marklines(1) marklines(1)],[0 100],'--','Color',[0.8 0.8 0.8])
    plot([marklines(2) marklines(2)],[0 100],'--','Color',[0.8 0.8 0.8])

    thisyear=whichyear(ii);
    for i=2:length(sim.p.nameGroup)
        % create grey area
        x_vector = [p.m(p.ixStart(i):p.ixEnd(i)), fliplr(p.m(p.ixStart(i):p.ixEnd(i)))];
        inBetween = [squeeze(Bosmo_photo_min(round(thisyear),:,i)), fliplr(squeeze(Bosmo_photo_max(round(thisyear),:,i)))];
        fill(x_vector', inBetween,'g','LineStyle','none','FaceColor',cm(i-1,:),'FaceAlpha',0.2);
        Bio=squeeze(Bosmo_photo_mean(round(thisyear),:,i))';
        %remove the points with no biomass so it does not clutter the plot
        minlim=0.5;
        [Bio]=removeTooSmallPoints(Bio,minlim);
        plot(p.m(p.ixStart(i):p.ixEnd(i)),Bio,'Color',cm(i-1,:),'LineWidth',1.5,'DisplayName',['Generalist \alpha=',num2str(alphaF_these(i-1))])
        
    end
    % wrap up axis and titel
    set(gca,'XScale','log');
    ax=gca;
    ax.Layer = 'top';
    ax.LineWidth=1;
    ax.FontSize=12;
    xlim([min(p.mLower(p.ixStart(2):p.ixEnd(2))) max(p.mUpper(p.ixStart(2):p.ixEnd(2)))])
    ylim([0 100])
    ax.YAxis.TickLabels={'0';'';'100'};
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