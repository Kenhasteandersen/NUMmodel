% It loads NPP matrices for different parameters at different water columns
% and plots barplots of NPP under varying parameters

% All parameters vectors
newRemin2= logspace(log10(0.01), log10(0.5), 10);
newMortHTL = logspace(log10(0.05), log10(0.9), 10);
newReminPOM= logspace(log10(0.01), log10(0.5), 10);
newSinkingPOM = logspace(log10(0.05), log10(0.9), 10);
newFracHTL= logspace(log10(0.1), log10(1), 10);

params={'remin2','mortHTL','reminPOM','sinkingPOM','fracHTL_to_N'};
directory_nm='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Sensitivity analysis\';
lat=[0;-5;24;55];

ilat=1;   % idx of selected latitude from 'lat'
iparam=2; % idx of selected parameter from 'params'


lat_to_find=lat(ilat);
param_vectors=[newRemin2;newMortHTL; newReminPOM; newSinkingPOM; newFracHTL];
param_val=param_vectors(iparam,end);

%% Figures for NPP - load saved matrices
ilat=1;
figure(ilat)
clf(ilat)
set(gcf, 'Color', 'white')
tiledlayout(length(lat),length(params));
for ilat=1:length(lat)
    for iparam=1:length(params)
        file_name=['NPP_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_nm,file_name]);
        NPP_annual_mean=mean(monthly_NPP_mean,2);
    nexttile
        boxplot(monthly_NPP_mean')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('NPP (mgC/m^{-2}/day)')
        ylim([0,3000])
        ttitle=['NPP at lat:',num2str(lat(ilat)),' for varying ',params{iparam}];
        title(ttitle)
    end
end

%% Figures for zremin
% loads saved matrices of annual median remineralization depth at 4 water columns
% for different parameters
directory='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\';

figure(ilat+1)
clf(ilat+1)
set(gcf, 'Color', 'white')
tiledlayout(length(lat),length(params));
for ilat=1:length(lat)
    for iparam=1:length(params)
        file_name=['zremin_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory,file_name]);
        zremin_annual_mean=mean(zremin,2);
    nexttile
        boxplot(zremin')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('NPP (mgC/m^{-2}/day)')
        ttitle=['z_{remin} at lat:',num2str(lat(ilat)),' for varying ',params{iparam}];
        title(ttitle)
    end
end
%%   Still zremin
figure(ilat+2)
clf(ilat+2)
set(gcf, 'Color', 'white')
tiledlayout(length(lat),length(params));
for ilat=1:length(lat)
    for iparam=1:length(params)
        file_name=['zremin_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory,file_name]);
        zremin_annual_mean=mean(zremin,2);
    nexttile
        plot(-zremin_annual_mean')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('NPP (mgC/m^{-2}/day)')
        ttitle=['z_{remin} at lat:',num2str(lat(ilat)),' for varying ',params{iparam}];
        title(ttitle)
        axis tight
    end
end
%% Figures for Bphyto
% loads saved matrices of annual median phytoplankton-0.6 biomass
% at 4 water columns
% for different parameters
directory='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\Sensitivity analysis\';

figure(ilat+3)
clf(ilat+3)
set(gcf, 'Color', 'white')
tiledlayout(length(lat),length(params));
for ilat=1:length(lat)
    for iparam=1:length(params)
        file_name=['Bphyto_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory,file_name]);
        Bphyto_annual_mean=mean(Bph_param,2);
    nexttile
        boxplot(Bph_param')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('B_{phyto} (\mugC/l)')
        ttitle=['B_{phyto0.6} at lat:',num2str(lat(ilat)),' for varying ',params{iparam}];
        title(ttitle)
    end
end
%% Save Btot, which will be used to calculate monthly averages

% load sim in for loop
% Compile Btot for each param
Btot=zeros(length(newReminPOM),25*365,23);
iparam=2;
for ixparam=1:length(newReminPOM)
    file_name=['sim',params{iparam},'_',num2str(param_vectors(iparam,ixparam)),...
        'lat',num2str(lat(ilat)),'.mat'];
    load([directory_nm,file_name]);
    Btot(ixparam,:,:)=squeeze(sum(sim.B(:,:,:),3));
end

save(['Btot',params{iparam},'_lat',num2str(lat_to_find),'.mat'],'Btot');    
%% Saving matrices of averages Ntot, Btot, zremin_monthly_mean
noYears=25;
depth_layers=[23;21;21;19]; %maunally identified
%
for ilat=1:length(lat)
monthly_Btot=zeros(length(newRemin2),noYears,depth_layers(ilat),12);
monthly_Btot_mean=zeros(length(newRemin2),depth_layers(ilat),12);
monthly_Ntot=monthly_Btot;
monthly_Ntot_mean=monthly_Btot_mean;
Btot=zeros(length(newRemin2),noYears*365,depth_layers(ilat));
Ntot=Btot;
for iparam=1:length(params)              % iterate for parameter
    for ixparam = 1:length(newRemin2)    % iterate for different values of the specific parameter
        file_name=['sim',params{iparam},'_',num2str(param_vectors(iparam,ixparam)),...
        'lat',num2str(lat(ilat)),'.mat'];
        load([directory_nm,file_name]);
        Btot(ixparam,:,:)=squeeze(sum(sim.B(:,:,:),3));
        Ntot(ixparam,:,:)=squeeze(sim.N(:,:));
        for i=noYears-3:noYears
            for idepth=1:depth_layers(ilat)
                monthly_Btot(ixparam,i,idepth,:)=reshapeCellToArrayAvg(squeeze(Btot(ixparam,:,idepth)),i);
                monthly_Ntot(ixparam,i,idepth,:)=reshapeCellToArrayAvg(squeeze(Ntot(ixparam,:,idepth)),i);
            end
        end
  % monthly Btot averaged over the last 3 years of the simulation
    monthly_Btot_mean(iparam,:,:)=squeeze(mean(monthly_Btot(ixparam,noYears-3:noYears,:,:),2));
    monthly_Ntot_mean(iparam,:,:)=squeeze(mean(monthly_Ntot(ixparam,noYears-3:noYears,:,:),2));
    end
save(['Btot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Btot_mean');    
save(['Ntot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Ntot_mean');    
end
end
%
for ilat=1:length(lat)
    for iparam=1:length(params) 
        file_name=['Ntot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load(file_name);
        for ixparam=1:length(newReminPOM)
            for j=1:12
                [~,id(ixparam,:,j)]=max(abs(diff(monthly_Ntot_mean(ixparam,:,j))));
                zremin(ixparam,j)=sim.z(id(ixparam,:,j)+1);
            end
        end
        save(['zremin_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'zremin');    
    end
end
%% PLOTS of vertical profiles for Ntot & Btot
vr=squeeze(monthly_Ntot_mean(i,:,j));
% [~, id]=
zremin=sim.z(id+1);

figure(3)
clf(3)
set(gcf,'Color','white')
tiledlayout(length(lat),12);

for i=2:5%1:length(newReminPOM)
    for j=1:12
nexttile
    plot(squeeze(monthly_Ntot_mean(i,:,j)),-sim.z,'LineWidth',2)
    hold on 
    plot(squeeze(monthly_Btot_mean(i,:,j)),-sim.z,'LineWidth',2)
    xlabel('concentration')
    ylabel('depth')
    % legend('N','Btot') 
    end
end
