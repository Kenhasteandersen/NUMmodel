%
% Sensitivity analysis of a Water column simulation for specific
% coordinates.
% First extract NPP from satellite data by running
% EXTRACT_files_NPP_latlon.m
%
% load('NPP_extracted_seasonal.mat')
% lat_to_find=55; lon_to_find=-40;name_png='sensitivityNPP';noYears=5;

function Sensitivity_mortHTL(lat_to_find,lon_to_find,name_png,NPP_extracted,noYears)
   
saved_png = append(name_png,'.png');

mHTL = 1;
mortHTL = .15;
bHTLdecline = true;
bHTLquadratic = true;   
runningTime = noYears*365;
sinkingPOM = 0.78622;
newMortHTL = logspace(log10(0.05), log10(0.9), 10);
newParameter=newMortHTL;
param_str='mortHTL';
param_str4plots=param_str;


Ntotal = zeros(length(newParameter), runningTime);
NPP = zeros(length(newParameter), runningTime);
NPP_annual = zeros(length(newParameter));

for  iparam = 1:length(newParameter)
     
        p = parametersWatercolumn( setupNUMmodel );
        setHTL(newMortHTL(iparam), mHTL, bHTLquadratic, bHTLdecline)
        setSinkingPOM(p, sinkingPOM)
        p.tEnd = runningTime;
        sim = simulateWatercolumn(p,lat_to_find,lon_to_find);
        simWCSeason = calcFunctionsTh(sim);
        save(['sim',param_str,'_',num2str(newParameter(iparam)),'lat',num2str(lat_to_find),'.mat'],'sim');    
        save(['simWC',param_str,'_',num2str(newParameter(iparam)),'lat',num2str(lat_to_find),'.mat'],'simWCSeason');    
        Ntotal(iparam,:) = simWCSeason.Ntot;
        NPP(iparam,:) = simWCSeason.ProdNet; % (mgC / m2 / day)
        Bph_orig(iparam,:,:,:) = simWCSeason.Bph;
        % NPP_annual(i) = simWCSeason.ProdNetAnnual;
end

%%
final_years=3;
monthly_NPP=zeros(length(newParameter),noYears,12);
monthly_NPP_mean=zeros(length(newParameter),12);
monthly_Bphyto=zeros(length(newParameter),noYears,12,length(sim.z));
monthly_Bph_mean=zeros(length(newParameter),12,length(sim.z));

for iparam = 1:length(newParameter)
    Bph_allSizes=squeeze(sum(Bph_orig(iparam,:,:,:),4));
    for i=noYears-3:noYears
         monthly_NPP(iparam,i,:)=reshapeCellToArrayAvg(NPP(iparam,:),i);
         for k=1:length(sim.z)
            monthly_Bphyto(iparam,i,:,k)=reshapeCellToArrayAvg(Bph_allSizes(:,k),i);
         end
    end
  % monthly NPP averaged over the last 3 years of the simulation
    monthly_NPP_mean(iparam,:)=squeeze(mean(monthly_NPP(iparam,noYears-3:noYears,:),2));
    monthly_Bph_mean(iparam,:,:)=squeeze(mean(monthly_Bphyto(iparam,noYears-3+1:noYears,:),2));
end
% Save matrix with monthly NPP averaged over the last 3 years, for every
% parameter
save(['NPP_monthly_mean_',param_str,'lat',num2str(lat_to_find),'.mat'],'monthly_NPP_mean');    
save(['Bph_monthly_mean_',param_str,'lat',num2str(lat_to_find),'.mat'],'monthly_Bph_mean');    

%% -------------------------------------------
%                  PLOTS
%-------------------------------------------

figure(1)
clf(1)
tiledlayout(1,2)
nexttile
for i = 1:length(newParameter)
        plot(monthly_NPP_mean(i,:), 'o-r', 'LineWidth', 1)
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

legend('NUM_{last3year.avg}','Eppley model', 'Standard VGPM', 'CAFE','Location','best')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
my_title = append(param_str4plots,' = ', string(newParameter(i)));
mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find),', ', my_title);
title(mTitle)
end

nexttile
geoscatter(lat_to_find, lon_to_find, 'filled')
exportgraphics(gcf,[saved_png])       
  
saved_png2 = append(name_png,'_details.png');
                
figure(2)
clf(2)
set(gcf,'Position',get(0,'ScreenSize'))
t=tiledlayout(2,round(length(newParameter)/2));
for iparam = 1:length(newParameter)

    nexttile                                   
        plot(squeeze(monthly_NPP_mean(iparam,:)), 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')
    xlabel('Time (month)')
    ylabel('NPP (mgC / m^2 /day)')
    mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
    my_title = append(param_str4plots,' = ', string(newParameter(iparam)));
    title(my_title)
end
title(t,mTitle, 'FontSize', 24)
lg  = legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','NumColumns',2); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

exportgraphics(gcf,[saved_png2])       

end