% load matrix of gloabl simulation
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\simGlobal.mat')
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\simGlobalF.mat')
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Develop-matrices\simGlobal5yearsF.mat')
% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\simNUMFun.mat')
%%
% sim=simGlobalF;
% sim=simGlobal5yearsF;
% Total global biomass
%

figure(2)
clf
tiledlayout(3,1)

nexttile
cbar = panelGlobal(sim.x,sim.y, squeeze(log10(mean(sim.Bpico,1)/1000)), [-1, 1], sProjection='Eckert4');
cbar.Label.String = "g_C/m^2";
title('pico')

nexttile
cbar = panelGlobal(sim.x,sim.y, squeeze(log10(mean(sim.Bnano,1)/1000)), [-1, 1], sProjection='Eckert4');
cbar.Label.String = "g_C/m^2";
title('nano')

nexttile
cbar = panelGlobal(sim.x,sim.y, squeeze(log10(mean(sim.Bmicro,1)/1000)), [-1, 1], sProjection='Eckert4');
cbar.Label.String = "g_C/m^2";
title('micro')

%%
%*************************************************
%-------------------------------------------------
%            Annual NPP, Btotal, Chl-a
%-------------------------------------------------
%*************************************************

plotGlobalFunctions(sim);

%% -----------------------------------------------
%*************************************************
%-------------------------------------------------
%                      AMT
%-------------------------------------------------
%*************************************************

load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Model EVALUATION\matrices2022\dataprotistsbiomass_corrected.mat');
% biomass in gC/m^2
load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\simGlobalF.mat')
% biomass in (mgC/m2)

% sim=simGlobalF;
T   = dataprotistsbiomass; % import table with AMT data
x   = T(:,2); % longitude
y   = T(:,5); % latitude
B   = T(:,3); % Biomass of pico- and nano-plankton
transect  = table2array(T(:,4));
month_corrected=table2array(T(:,8));
% converting tables to arrays
lat_dat = table2array(x); 
lon_dat = table2array(y);
B_dat   = table2array(B);

% indices for the different shapes/months in Serra-Pompei et al. (2022)
romb     =  1:27;   % diamond/ AMT-12
romb_may =  1:14;    % becaus first half of that cruise was in May  
cuadr    =  28:53;  % square/ AMT-14
tri      =  54:77;  % triangle/ AMT-13
tri_sep  = tri(10:24);
% AMT-13 started from North to South, so we need to flip indexes
% tri      = flip(tri);
% tri_sep  = flip(tri_sep);
lonWrapped = wrapTo180(sim.x); % convert longitude ordinates from [0,360[-->[-180,180]

% Global simulation from NUMmodel 
% for may-june-october
yrNo = sim.p.tEnd/365;
may  = (yrNo-1)*12+5;
june = may+1;
oct  = (yrNo-1)*12+10;
sep  = oct-1;

%  
Bpico_5  =  sim.Bpico(may,:,:)/1000; % convert to gC/m^2 for all
Bpico_6  =  sim.Bpico(june,:,:)/1000;
Bpico_9 =  sim.Bpico(sep,:,:)/1000;
Bpico_10 =  sim.Bpico(oct,:,:)/1000;
Bnano_5  =  sim.Bnano(may,:,:)/1000;
Bnano_6  =  sim.Bnano(june,:,:)/1000;
Bnano_9 =  sim.Bnano(sep,:,:)/1000;
Bnano_10 =  sim.Bnano(oct,:,:)/1000;

% Taking monthly average biomass of the model data
Bg_5avg  =  squeeze(sum(Bpico_5+Bnano_5,1)/length(may));
Bg_6avg  =  squeeze(sum(Bpico_6+Bnano_6,1)/length(june));
Bg_9avg =  squeeze(sum(Bpico_9+Bnano_9,1)/length(sep));
Bg_10avg =  squeeze(sum(Bpico_10+Bnano_10,1)/length(oct));

%----------------------------------------------------------
% Find model Biomass for the same coordinates with the data
%----------------------------------------------------------
Pbiom_model5=zeros(1,length(lat_dat));
Pbiom_model6=zeros(1,length(lat_dat));
Pbiom_model9=zeros(1,length(lat_dat));
Pbiom_model10=zeros(1,length(lat_dat));


for i=1:length(lat_dat)

     [ ~, idx_lonG] = min( abs( lonWrapped-lon_dat(i) ) );

     [ ~, idx_latG] = min( abs( sim.y-lat_dat(i) ) );

     ix_lon(i)=idx_lonG;
     ix_lat(i)=idx_latG;

% next step is to find combined coordinates
    Pbiom_model6(i)=Bg_6avg(idx_lonG(1),idx_latG(1));
    Pbiom_model5(i)=Bg_5avg(idx_lonG(1),idx_latG(1));
    Pbiom_model10(i)=Bg_10avg(idx_lonG(1),idx_latG(1));
    Pbiom_model9(i)=Bg_9avg(idx_lonG(1),idx_latG(1));

end
Pbiom_modelAMT12=[Pbiom_model5(romb_may), Pbiom_model6((length(romb_may))+1:length(romb))];
Pbiom_modelAMT13=[Pbiom_model10((length(lat_dat)-length(tri)+1):62),Pbiom_model9(tri_sep) ];


% plot model vs data for each month
% cmap=[23, 25, 82; 64, 135, 163; 115, 209, 184]./255;
cmap=[169,223,205; 149,167,255; 135,73,163]./255; % mermaiddoc

figure(3)
clf(3)
subplot(1,4,2)
% axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
%     'Frame','on','Grid','on','MeridianLabel','off','ParallelLabel','on','fontsize',8)
hold on
plot(B_dat(romb),lat_dat(romb),'d','markerfacecolor',cmap(1,:));%,'filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
% hold on
% Pbiom_modelAMT12=[Pbiom_model5(romb_may), Pbiom_model6((length(romb_may))+1:length(romb))];
% plot(sim.y(idx_lat),Bg_6avg(idx_lon,idx_lat))
% plot(sim.y(ix_lat),Bg_6avg(ix_lon,ix_lat))
plot(Pbiom_model6(romb),sim.y(ix_lat(romb)),'LineWidth',2,'color',cmap(1,:))
% plot(sim.y(ix_lat(romb)),Pbiom_modelAMT12,'LineWidth',2)
hold off

ylabel('latitude')
xlabel("biomass (gC/m^2)")
title("AMT 12:May-June")

subplot(1,4,3)
plot(B_dat(tri),lat_dat(tri),"^",'markerfacecolor',cmap(2,:))
hold on
plot(Pbiom_modelAMT13,sim.y(ix_lat(tri)),'LineWidth',2,'color',cmap(2,:))
% ylabel('latitude')
xlabel("biomass (gC/m^2)")
title("AMT 13:Oct-Sep (Southward)")

subplot(1,4,4)
plot(B_dat(cuadr),lat_dat(cuadr),"s",'markerfacecolor',cmap(3,:))
hold on
plot(Pbiom_model5(cuadr),sim.y(ix_lat(cuadr)),'LineWidth',2,'color',cmap(3,:))
% ylabel('latitude')
xlabel("biomass (gC/m^2)")
title("AMT 14:May")


%---------------------
% Plot transects
%---------------------
% cmap=[23, 25, 82; 64, 135, 163; 115, 209, 184]./255;
latlim=[-47,60];
lonlim=[-48,-1];

subplot(1,4,1)
axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','fontsize',8)
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
SC1=scatterm(lat_dat(transect==1),lon_dat(transect==1),15,'d','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
SC2=scatterm(lat_dat(transect==3),lon_dat(transect==3),15,'^','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');
SC3=scatterm(lat_dat(transect==2),lon_dat(transect==2),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
SC3=scatterm(sim.y(ix_lat(54:77)),sim.x(ix_lon(54:77)),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');

leg=legend([SC1,SC2,SC3],{'AMT 12','AMT 13','AMT 14'},'fontsize',8);

leg.ItemTokenSize(1) = 10;
leg.Location='southoutside';
leg.NumColumns=3;
%%
%*************************************************
%-------------------------------------------------
%        COPEPOD BIOMASS - Lopez anadon
%-------------------------------------------------
%*************************************************

dat=importdata('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\data_lopez_anadon.xlsx');
sim=simGlobalF;
% integrated over depth - top 5 layers (top 60m)
simInt=calcIntegrateGlobal(simGlobalF,simGlobalF.B , false);
simIntegrate=simInt; % in mgC/m^2

y=sim.y; x=sim.x;

lat_dat=dat.data(:,10);
lon_dat=dat.data(:,9);
lon_dat=mod(lon_dat(:),360);
cop_biom_dat=dat.data(:,7);
naup_biom_dat=dat.data(:,8);
station_dat=dat.data(:,2);
cop_sfrac_dat=dat.data(:,1);
cop_body_mass=dat.data(:,5).*1e-3;
% conversion of body mass to length
cop_mm=(cop_body_mass./exp(-16.41)).^(1/2.74);%*1e-3;
cop_mm_mat=reshape(cop_mm,[length(cop_mm)/4,4]);
cop_biom_mat=reshape(cop_biom_dat,[length(cop_mm)/4,4]);
% wght_mean=sum((cop_biom_mat.*cop_mm_mat),2)./sum(cop_biom_mat,2);

s1=find(cop_sfrac_dat==1); %Smallest
s2=find(cop_sfrac_dat==2);
s3=find(cop_sfrac_dat==3);
s4=find(cop_sfrac_dat==4); %Largest

% find Copepod () size groups in the model
cop_group_ix1=find(sim.p.typeGroups==10);
cop_group_ix2=find(sim.p.typeGroups==11);
cop_group_ix = [cop_group_ix1 cop_group_ix2];

C_ix=sim.p.ixStart(cop_group_ix(1)):sim.p.ixEnd(cop_group_ix(end)); % indices of Copepod groups
sb4=find(sim.p.m(C_ix)>50);
sb3=find(sim.p.m(C_ix)>10 & sim.p.m(C_ix)<=50);
sb2=find(sim.p.m(C_ix)>1 & sim.p.m(C_ix)<=10);
sb1=find(sim.p.m(C_ix)<=1);
Cdeach=simIntegrate;% depth-integrated
% Cdeach= cell2mat(struct2cell(simIntegrate)); % this converts struct array to matrix

Cdtot_4D_1=sum(Cdeach(:,:,:,sb1),4); % adds all the small ones
Cdtot_4D_2=sum(Cdeach(:,:,:,sb2),4);
Cdtot_4D_3=sum(Cdeach(:,:,:,sb3),4);
Cdtot_4D_4=sum(Cdeach(:,:,:,sb4),4);


Cbiom_model=zeros(1,length(lat_dat(s1)));
Cbiom_model_1=zeros(1,length(lat_dat(s1)));
Cbiom_model_2=zeros(1,length(lat_dat(s1)));
Cbiom_model_3=zeros(1,length(lat_dat(s1)));
Cbiom_model_4=zeros(1,length(lat_dat(s1)));

yrNo = sim.p.tEnd/365;
oct  = (yrNo-1)*12+10;
sep  = oct-1;


for i=1:length(lat_dat(s1))
    idx1=find(y>lat_dat(i)-2.8 & y<lat_dat(i)+2.8);
    idx2=find(x>lon_dat(i)-3 & x<lon_dat(i)+3);
%     Cbiom_model(i)=nansum(mean(Cdtot_4D(idx2(1),idx1(1),:,9*30+14:10*30+10).*permute(deltaz,[3,2,1]),4),3);
% take the average for September-October
    Cbiom_model_1(i)=nansum(mean(Cdtot_4D_1(sep:oct,idx2(1),idx1(1)),1),1);
    Cbiom_model_2(i)=nansum(mean(Cdtot_4D_2(sep:oct,idx2(1),idx1(1)),1),1);
    Cbiom_model_3(i)=nansum(mean(Cdtot_4D_3(sep:oct,idx2(1),idx1(1)),1),1);
    Cbiom_model_4(i)=nansum(mean(Cdtot_4D_4(sep:oct,idx2(1),idx1(1)),1),1);
end

cop_tot_dat=cop_biom_dat(s1)+cop_biom_dat(s2)+cop_biom_dat(s3)+cop_biom_dat(s4);

% Plot comparing data and model

size_mat=[cop_biom_dat(s4),cop_biom_dat(s3),cop_biom_dat(s2),cop_biom_dat(s1)+naup_biom_dat(s1)]; %mgC/m^2
size_mat_mod=[Cbiom_model_4',Cbiom_model_3',Cbiom_model_2',Cbiom_model_1'];

x0=0;
y0=0;
width=18;
height=7;


fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
tiledlayout(1,3,'Padding','Compact');

nexttile(2)
scatter(0,0,'markeredgecolor','none')
hold on
hb=bar(station_dat(s1),size_mat_mod./1000,'stacked');
cmap=flip(brewermap(4,'BuPu'));

hb(1).FaceColor = cmap(1,:); %large
hb(2).FaceColor = cmap(2,:);
hb(3).FaceColor = cmap(3,:);
hb(4).FaceColor = cmap(4,:); %small

leg1=legend([hb(4) hb(3) hb(2) hb(1)],{'<1','1 - 10','10 - 50','50<'});
legend boxoff
title(leg1,'Body-mass range')
set(gcf,'color','w');
ylim([0 10])
title('  {\bfa.} Model','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Biomass [gC m^{-2}]')
set(gca, 'YDir','reverse')
ylabel('Station')

nexttile(3)
scatter(0,0,'markeredgecolor','none')
hold on
hb=barh(station_dat(s1),size_mat./1000,'stacked');
cmap=flip(brewermap(4,'BuPu'));

hb(1).FaceColor = cmap(1,:);
hb(2).FaceColor = cmap(2,:);
hb(3).FaceColor = cmap(3,:);
hb(4).FaceColor = cmap(4,:);
hold on
hp=plot(station_dat(s1),(Cbiom_model)./1000,'k--','linewidth',1);
leg2=legend([hb(4) hb(3) hb(2) hb(1)],{'0.4','3.3','13.6','86.8'});
legend boxoff
title(leg2,'Average body-mass')
xlim([0 10])
title('  {\bfb.} Data','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Biomass [gC m^{-2}]')
ylabel('Station')
set(gca, 'YDir','reverse')

c_p= [236 197 68]./255;
latlim = [-50 60];
lonlim = [-50 10];
nexttile(1)

axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on')
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
scatterm(lat_dat,lon_dat,10,'filled','markerfacecolor',c_p,'markeredgecolor','k')
text(lat_dat+2,lon_dat+2, num2str(1),'fontsize',8);
for i=1:24
    if i<10
        textm(lat_dat(i)+1,lon_dat(i)-5, num2str(i),'fontsize',6);
    else
        textm(lat_dat(i)+1,lon_dat(i)+2, num2str(i),'fontsize',6);        
    end

end

pos = get(gca, 'Position');
pos(1) = -0.45;
pos(3)=1.3;
set(gca, 'Position', pos)
set(gcf,'color','w');
