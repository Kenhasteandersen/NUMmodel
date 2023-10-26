% It loads NPP matrices for different parameters at different water columns
% and plots barplots of NPP under varying parameters
cd 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab'

% All parameters vectors
newSinkingPOM = logspace(log10(1), log10(50), 10);
params={'sinkingPOM'};
directory_nm='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\sinkingPOM sensitivity\';
directory_zremin='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\';
directory_h='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\';

lat=[0;-5;24;55];

% lat_to_find=lat(ilat);
param_vectors=newSinkingPOM;
noYears=25;
depth_layers=[23;21;21;19]; %maunally identified 
iparam=1;
%% CALCULATE zremin
% %Saving matrices of averages Ntot, Btot, zremin_monthly_mean
z={};
for ilat=1:length(lat)
% Initialize matrices
%                         param            years      depth         months
 monthly_Btot=zeros(length(newSinkingPOM),noYears,depth_layers(ilat),12);
 monthly_Ntot=monthly_Btot; monthly_Bcop=monthly_Btot; monthly_Buni=monthly_Btot;
 monthly_Bph=zeros(length(newSinkingPOM),noYears,depth_layers(ilat),12);
%                          param                   depth          months 
 monthly_Btot_mean=zeros(length(newSinkingPOM),depth_layers(ilat),12);
 monthly_Ntot_mean=monthly_Btot_mean;  monthly_Bcop_mean=monthly_Btot_mean; monthly_Buni_mean=monthly_Btot_mean;
 monthly_Bph_mean=zeros(length(newSinkingPOM),depth_layers(ilat),12);
 
 Btot=zeros(length(newSinkingPOM),noYears*365,depth_layers(ilat));
 Ntot=Btot;  Bcop=Btot;  Buni=Btot; 

    for ixparam = 1:length(newSinkingPOM)    % iterate for different values of the specific parameter
    % ixparam=1;
        file_name=['simWC',params{iparam},'_',num2str(param_vectors(iparam,ixparam)),...
        'lat',num2str(lat(ilat)),'.mat'];
        load([directory_nm,file_name]);
        z{ilat}=simWCSeason.z;
       %.............. find indices of Multicellular Plankton .............
        [~,idxM]=find((simWCSeason.p.typeGroups>=10) & (simWCSeason.p.typeGroups<100));
        ixMPlankton=(simWCSeason.p.ixStart(idxM(1)):simWCSeason.p.ixEnd(idxM(end)))-simWCSeason.p.idxB+1;
        [~,idxU]=find((simWCSeason.p.typeGroups<10));
        ixUPlankton=(simWCSeason.p.ixStart(idxU(1)):simWCSeason.p.ixEnd(idxU(end)))-simWCSeason.p.idxB+1;
       %...................................................................
        %     param, t,z                            t,z,size
        Btot(ixparam,:,:)=squeeze(sum(simWCSeason.B(:,:,1:end-1),3));
        Bcop(ixparam,:,:)=squeeze(sum(simWCSeason.B(:,:,ixMPlankton),3));
        Buni(ixparam,:,:)=squeeze(sum(simWCSeason.B(:,:,ixUPlankton),3));
        %                                       t,z 
        Ntot(ixparam,:,:)=squeeze(simWCSeason.N(:,:));
        Bph_all=sum(simWCSeason.Bph,3);

        for iyear=noYears-3:noYears
            for idepth=1:depth_layers(ilat)    %12 months
                monthly_Btot(ixparam,iyear,idepth,:)=reshapeCellToArrayAvg(squeeze(Btot(ixparam,:,idepth)),iyear);
                monthly_Bcop(ixparam,iyear,idepth,:)=reshapeCellToArrayAvg(squeeze(Bcop(ixparam,:,idepth)),iyear);
                monthly_Buni(ixparam,iyear,idepth,:)=reshapeCellToArrayAvg(squeeze(Buni(ixparam,:,idepth)),iyear);
                monthly_Bph(ixparam,iyear,idepth,:)=reshapeCellToArrayAvg(squeeze(Bph_all(:,idepth)),iyear);
                monthly_Ntot(ixparam,iyear,idepth,:)=reshapeCellToArrayAvg(squeeze(Ntot(ixparam,:,idepth)),iyear);
            end
        end
  % monthly Btot averaged over the last 4 years of the simulation
    monthly_Btot_mean(ixparam,:,:)=squeeze(mean(monthly_Btot(ixparam,noYears-3:noYears,:,:),2));
    monthly_Bcop_mean(ixparam,:,:)=squeeze(mean(monthly_Bcop(ixparam,noYears-3:noYears,:,:),2));
    monthly_Buni_mean(ixparam,:,:)=squeeze(mean(monthly_Buni(ixparam,noYears-3:noYears,:,:),2));
    monthly_Bph_mean(ixparam,:,:)=squeeze(mean(monthly_Bph(ixparam,noYears-3:noYears,:,:),2));
    monthly_Ntot_mean(ixparam,:,:)=squeeze(mean(monthly_Ntot(ixparam,noYears-3:noYears,:,:),2));
        for j=1:12
            [~,id(ixparam,:,j)]=max((monthly_Ntot_mean(ixparam,:,j)));
            [~,id_diff(ixparam,:,j)]=max(abs(diff(monthly_Ntot_mean(ixparam,:,j))));
            id_lat{ilat}=id;
            id_lat_diff{ilat}=id_diff;
            zreminPOM(ilat,ixparam,j)=z{ilat}(id_lat_diff{ilat}(ixparam,j)+1);
        end
       
    end
  save(['Btot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Btot_mean'); 
  save(['Bcop_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Bcop_mean');    
  save(['Buni_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Buni_mean');    
  save(['Bph_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Bph_mean');    
  save(['Ntot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'],'monthly_Ntot_mean');    
  save(['zremin_monthly_mean.mat'],'zreminPOM');    

end

%% FIGURES OF NPP, Bphyto & zremin for varying sinkingPOM
figure(1)
clf(1)
set(gcf, 'Color', 'white')
tiledlayout(length(lat),4);
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
        ttitle=['NPP at lat:',num2str(lat(ilat))];
        title(ttitle)

        file_name=['Bph_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_nm,file_name]);
        Bph_annual_mean_depthInt=squeeze(mean(monthly_Bph_mean,3));
    nexttile
        boxplot(Bph_annual_mean_depthInt')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('B_{phyto} (\mugC/l)')
        % ylim([0,3000])
        ttitle=['B_{phyto0.6} depth int at lat:',num2str(lat(ilat))];
        title(ttitle)
        
        % file_name=['zremin_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        % load([directory_zremin,file_name]);
        % zremin_annual_mean=mean(squeeze(zremin),2);
    nexttile
        % plot(-zremin_annual_mean')
        plot(-mean(zreminPOM(ilat,:,:),3)') %
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('Depth (m)')
        ttitle=['z_{remin} annual avg at lat:',num2str(lat(ilat))];
        title(ttitle)
        axis tight

     % file_nameNPP=['NPP_extracted',num2str(lat(ilat)),'lat.mat'];
     % load([directory_h,file_nameNPP]);

    file_nameNPPcafe=['NPP_extracted_CAFElat',num2str(lat(ilat))];
    load([directory_h,file_nameNPPcafe]);


nexttile                              % without this
        plot(monthly_NPP_mean(1,:), 'o-r', 'LineWidth', 1)% mgC/m2/day
        % plot(NPP_cell_month_mean{i, j}, 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
        plot(monthly_NPP_mean(2,:), 's-c', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(3,:), 'd-y', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(4,:), 'o:r', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(5,:), 's:c', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(6,:), 'd:y', 'LineWidth', 2)% mgC/m2/day
        plot(NPP_extracted(1,:), 'ob--')
        % plot(NPP_extracted(2,:), 'og--')
        % plot(NPP_extracted(3,:), 'om--')

    NUM1=num2str(round(param_vectors(iparam,1),2));
    NUM2=num2str(round(param_vectors(iparam,2),2));
    NUM3=num2str(round(param_vectors(iparam,3),2));
    NUM4=num2str(round(param_vectors(iparam,4),2));
    NUM5=num2str(round(param_vectors(iparam,5),2));
    NUM6=num2str(round(param_vectors(iparam,6),2));

    xlabel('Time (month)')
    ylabel('NPP (mgC / m^2 /day)')
    mTitle = append('Lat: ', string(lat(ilat)));
    title(mTitle)
    
    end
end
% lg  = legend(['NUM',params{iparam},NUM1],['NUM',params{iparam},NUM2],['NUM',params{iparam},NUM3],...
%     ['NUM',params{iparam},NUM4],['NUM',params{iparam},NUM5],['NUM',params{iparam},NUM6],...
%     'Eppley model', 'Standard VGPM', 'CAFE','NumColumns',3); 
lg  = legend(['NUM',params{iparam},NUM1],['NUM',params{iparam},NUM2],['NUM',params{iparam},NUM3],...
    ['NUM',params{iparam},NUM4],['NUM',params{iparam},NUM5],['NUM',params{iparam},NUM6],...
     'CAFE','NumColumns',3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout
%%
figure(2)
clf(2)
set(gcf, 'Color', 'white')
tiledlayout(length(lat),4);
for ilat=1:length(lat)
    for iparam=1:length(params)
        file_name=['Bph_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_nm,file_name]);
        Bph_annual_mean_depthInt=squeeze(mean(monthly_Bph_mean,3));
    nexttile
        boxplot(Bph_annual_mean_depthInt')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('B_{phyto} (\mugC/l)')
        % ylim([0,3000])
        ttitle=['B_{phyto0.6} depth int at lat:',num2str(lat(ilat))];
        title(ttitle)


        file_name=['Bcop_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_h,file_name]);
   nexttile
        boxplot(squeeze(mean(monthly_Bcop_mean(:,:,:),2))')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('B_{cop} (\mugC/l)')
        % ylim([0,3000])
        ttitle=['B_{cop} depth int at lat:',num2str(lat(ilat))];
        title(ttitle)

        file_name=['Btot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_h,file_name]);
    nexttile
        boxplot(squeeze(mean(monthly_Btot_mean(:,:,:),2))')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('B_{tot} (\mugC/l)')
        ttitle=['B_{tot} depth int at lat:',num2str(lat(ilat))];
        title(ttitle)

        file_name=['Buni_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_h,file_name]);
   nexttile
        boxplot(squeeze(mean(monthly_Buni_mean(:,:,:),2))')
        xticks(1:10)
        xticklabels(string(round(param_vectors(iparam,:),2))) % rounded to the 2nd decimal
        xlabel(string(params{iparam}))
        ylabel('B_{uni} (\mugC/l)')
        ttitle=['B_{uni} depth int at lat:',num2str(lat(ilat))];
        title(ttitle)
            
    end
end

%%
figure(3)
clf(3)
% set(gcf, 'Position', get(0, 'Screensize'))
set(gca,'color','w');

t=tiledlayout(2,2);
for ilat = 1:length(lat)
    file_nameNPP=['NPP_extracted',num2str(lat(ilat)),'lat.mat'];
    load([directory_h,file_nameNPP]);
    file_name=['NPP_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
    load([directory_nm,file_name]);
nexttile                              % without this
        plot(monthly_NPP_mean(1,:), 'o-r', 'LineWidth', 1)% mgC/m2/day
        % plot(NPP_cell_month_mean{i, j}, 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
        plot(monthly_NPP_mean(2,:), 's-c', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(3,:), 'd-y', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(4,:), 'o:r', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(5,:), 's:c', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(6,:), 'd:y', 'LineWidth', 2)% mgC/m2/day

        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

ylim([0,3000])
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
mTitle = append('Lat: ', string(lat(ilat)));
title(mTitle)
NUM1=num2str(round(param_vectors(iparam,1),2));
NUM2=num2str(round(param_vectors(iparam,2),2));
NUM3=num2str(round(param_vectors(iparam,3),2));
NUM4=num2str(round(param_vectors(iparam,4),2));
NUM5=num2str(round(param_vectors(iparam,5),2));
NUM6=num2str(round(param_vectors(iparam,6),2));
end
lg  = legend(['NUM',params{iparam},NUM1],['NUM',params{iparam},NUM2],['NUM',params{iparam},NUM3],...
    ['NUM',params{iparam},NUM4],['NUM',params{iparam},NUM5],['NUM',params{iparam},NUM6],...
    'Eppley model', 'Standard VGPM', 'CAFE','NumColumns',3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

%%  
figure(4)
clf(4)
% set(gcf, 'Position', get(0, 'Screensize'))
set(gcf,'color','w');

t=tiledlayout(2,2);
for ilat = 1:length(lat)
file_nameNPPcafe=['NPP_extracted_CAFElat',num2str(lat(ilat))];
    load([directory_h,file_nameNPPcafe]);
    file_name=['NPP_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
    load([directory_nm,file_name]);
nexttile                              % without this
        plot(monthly_NPP_mean(1,:), 'o-r', 'LineWidth', 1)% mgC/m2/day
        % plot(NPP_cell_month_mean{i, j}, 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
        plot(monthly_NPP_mean(2,:), 's-c', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(3,:), 'd-y', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(4,:), 'o:r', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(5,:), 's:c', 'LineWidth', 2)% mgC/m2/day
        plot(monthly_NPP_mean(6,:), 'd:y', 'LineWidth', 2)% mgC/m2/day

        plot(NPP_extracted(1,:), 'ob--')

ylim([0,3000])
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
mTitle = append('Lat: ', string(lat(ilat)));
title(mTitle)
NUM1=num2str(round(param_vectors(iparam,1),2));
NUM2=num2str(round(param_vectors(iparam,2),2));
NUM3=num2str(round(param_vectors(iparam,3),2));
NUM4=num2str(round(param_vectors(iparam,4),2));
NUM5=num2str(round(param_vectors(iparam,5),2));
NUM6=num2str(round(param_vectors(iparam,6),2));
end
lg  = legend(['NUM',params{iparam},NUM1],['NUM',params{iparam},NUM2],['NUM',params{iparam},NUM3],...
    ['NUM',params{iparam},NUM4],['NUM',params{iparam},NUM5],['NUM',params{iparam},NUM6],...
 'CAFE','NumColumns',3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

%% PLOTS of vertical profiles for Ntot & Btot

figure(22)
clf(22)
set(gcf,'Color','white')
tiles=tiledlayout(4,10);
for ilat=1:length(lat)
    for ixparam=1:length(newSinkingPOM)
        % for j=1:12
        file_name1=['Ntot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_zremin,file_name1]);       
        file_name2=['Btot_monthly_mean_',params{iparam},'lat',num2str(lat(ilat)),'.mat'];
        load([directory_zremin,file_name2]);       

nexttile
        j=7; %months
        plot(squeeze(monthly_Ntot_mean(ixparam,:,j)),-z{ilat},'LineWidth',2)
        hold on 
        plot(squeeze(monthly_Btot_mean(ixparam,:,j)),-z{ilat},'LineWidth',2)
        plot(-z{ilat}(id_lat{ilat}(ixparam,j)),'gd')
        plot(-z{ilat}(id_lat_diff{ilat}(ixparam,j)),'ms')
        plot(-zreminPOM(ilat,ixparam,j),'yo')
        xlabel('concentration')
        ylabel('depth')
        ylim([-2000,0])
        % legend('N','Btot','z_{remin}') 
     end
    % end
    mTitle=['Monthly Vertical profiles at lat:',num2str(lat(ilat)),' for varying ',params{iparam},', month:',num2str(j)];
    title(tiles,mTitle, 'FontSize', 20)
end
