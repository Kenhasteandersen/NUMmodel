%
% Sensitivity analysis of a Water column simulation for specific
% coordinates.
% First extract NPP from satellite data by running
% EXTRACT_files_NPP_latlon.m
%
%%
%--------------------------------------------
%        Initialize coordinates
%--------......................--------
        lat_to_find = 24;
        lon_to_find = -158;
noYears=10;
name_png = 'sensitivityNPP_';
directory = 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\Comparisons with data\';
    
saved_png = append(directory,name_png,'.png');

mHTL = 1;
mortHTL = .15;
bHTLdecline = true;
bHTLquadratic = true;   
runningTime = noYears*365;
sinkingPOM = 30;
newStages = logspace(log10(0.05), log10(0.9), 10);
newPops = logspace(log10(0.05), log10(0.9), 10);




nPop = 4:10; nCopepods = 4:10;
n = 10; nPOM = 1;
% mAdultPassive = logspace(log10(0.2), log10(0.5), floor(nPop/2));
% mAdultActive = logspace(log10(1), log10(1000), floor(nPop/2) + mod(nPop,2));
% 
newParameter=nPop;
param_str='nStages';
param_str4plots=param_str;

%--------------------------------------------


for i =1:length(nPop)
   for j = 1:length(nCopepods) % life stages
    mAdultPassive = logspace(log10(0.2), log10(0.5), floor(nPop(i)/2));
    mAdultActive = logspace(log10(1), log10(1000), floor(nPop(i)/2) + mod(nPop(i),2));
        p = setupNUMmodel(mAdultPassive, mAdultActive, n,nCopepods(j),nPOM);
        p = parametersWatercolumn(p); 
        setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
        setSinkingPOM(p, sinkingPOM);
        p.tEnd = runningTime;
        sim = simulateWatercolumn(p,lat_to_find,lon_to_find);
        simWCSeason = calcFunctionsTh(sim);
        save(['sim',param_str,'_',num2str(i),'stage',num2str(j),'lat',num2str(lat_to_find),'.mat'],'sim');    
        save(['simWCnPop',num2str(i),'_nStages',num2str(j),'lat',num2str(lat_to_find),'.mat'],'simWCSeason');    
        Ntotal(i,j,:) = simWCSeason.Ntot;
        NPP(i,j,:) = simWCSeason.ProdNet; % (mgC / m2 / day)
        Bph_orig(i,j,:,:,:) = simWCSeason.Bph;
        % NPP_annual(i) = simWCSeason.ProdNetAnnual;    

    end
end

%%
final_years=3;
monthly_NPP=zeros(length(newParameter),length(nCopepods),noYears,12);
monthly_NPP_mean=zeros(length(newParameter),length(nCopepods),12);
monthly_Bphyto=zeros(length(newParameter),length(nCopepods),noYears,12,length(sim.z));
monthly_Bph_mean=zeros(length(newParameter),length(nCopepods),12,length(sim.z));

for iparam = 1:length(newParameter)
    for iparam2 = 1:length(nCopepods)
        Bph_allSizes=squeeze(sum(Bph_orig(iparam,iparam2,:,:,:),5));
        for i=noYears-3:noYears
         monthly_NPP(iparam,iparam2,i,:)=reshapeCellToArrayAvg(NPP(iparam,iparam2,:),i);
         for k=1:length(sim.z)
            monthly_Bphyto(iparam,iparam2,i,:,k)=reshapeCellToArrayAvg(Bph_allSizes(:,k),i);
         end
        end
  % monthly NPP averaged over the last 3 years of the simulation
    monthly_NPP_mean(iparam,iparam2,:)=squeeze(mean(monthly_NPP(iparam,iparam2,noYears-3:noYears,:),3));
    monthly_Bph_mean(iparam,iparam2,:,:)=squeeze(mean(monthly_Bphyto(iparam,iparam2,noYears-3+1:noYears,:,:),3));
    end
end
% Save matrix with monthly NPP averaged over the last 3 years, for every
% parameter
save(['NPP_monthly_mean_',param_str,'lat',num2str(lat_to_find),'.mat'],'monthly_NPP_mean');    
save(['Bph_monthly_mean_',param_str,'lat',num2str(lat_to_find),'.mat'],'monthly_Bph_mean');    

%%
% load('NPP_extracted.mat')
% load('NPP_extracted55lat.mat')

figure(1)
clf(1)
set(gcf,'color','w');

tiledlayout(1,3)

nexttile
for i = 1:length(nPop)

    for j = 1:length(nCopepods)

        plot(squeeze(monthly_NPP_mean(i,j,:)), 'o-r', 'LineWidth', 1)
        % plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)% mgC/m2/day
        % % plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','Location','best')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
 mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
 my_title = append('#pop = ', string(nPop(i)), ' / stages = ', string(nCopepods(j)));
        title(mTitle)

    end
end

nexttile
s=surface(nPop,nCopepods,squeeze(monthly_NPP_mean(:,:,4))) % in April, when it goes off 
colorbar
axis tight
% plotSizespectrum(simEutro);
s.EdgeColor = 'none';
xlabel('copepod populations')
ylabel('life stages')
title('mean NPP in April')


nexttile
s=surface(nPop,nCopepods,squeeze(mean(monthly_NPP_mean,3))) % in April, when it goes off 
colorbar
axis tight
% plotSizespectrum(simEutro);
s.EdgeColor = 'none';
xlabel('copepod populations')
ylabel('life stages')
title('annual mean NPP ')
colormap(whitejet);                                               
colorbar;
% geoscatter(lat_to_find, lon_to_find, 'filled')
% 
% save('NPP_cell_month_mean.mat');

% exportgraphics(gcf,['C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\...' ...
    % 'Trophic efficiency\Compare Biomass\Comparisons with data\NPPsensitivity1.png'])

% exportgraphics(gcf,[saved_png])       


%%
   
saved_png2 = append(directory,name_png,'_details.png');
% uplim=max(max(NPP_cell_month_mean),max(NPP_extracted));
                
figure(2)
clf(2)
t=tiledlayout(2,4)
for i = 1:length(nPop)

    for j = 1:length(nStages)

nexttile
        plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)% mgC/m2/day
        % plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

% legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','Location','best')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
 mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
 my_title = append('#Pop = ', string(nPop(i)), ' / nStages = ', string(nStages(j)));
        title(my_title)

    end
end
title(t,mTitle, 'FontSize', 24)
lg  = legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','NumColumns',2); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout


exportgraphics(gcf,[saved_png2])       

%%  GROUP CONTRIBUTIONS
depth_layer=1;
Bnum=squeeze(mean(sim.B(:,depth_layer,:),1)); % average Biomass at specific depth layer

PicoNanoMicroBarplots(sim,Bnum,depth_layer)

exportgraphics(gcf,[append(directory,name_png,'_layer',num2str(depth_layer),'_barplots.png')])       






