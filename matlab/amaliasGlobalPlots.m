        figure(1)
        clf
        plotGlobalAmalia(sim); 


        figure(2)
        clf
        plotGlobalFunctions(sim); 

        % lat = -5;
        % lon = 5;
        lat=[0;-5;24;55];
        lon=[172;5;-158;-40];
        
        for ilat=1:length(lat)
            figure
            plotWatercolumnTime(sim,lat(ilat),lon(ilat), depthMax=200);

            figure
            plotSizespectrumTime(sim,1,lat(ilat),lon(ilat));
            title(sprintf('Size spectrum at (%3.0f,%3.0f).\n',[lat(ilat),lon(ilat)]))
        end

%%
directory_nm='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\sinkingPOM sensitivity\';
directory_h='C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Documents\GitHub\NUMmodel\matlab\';

lat=[0;-5;24;55];
lon=[172;5;-158;-40];
% lat_to_find=lat(ilat);
noYears=5;
depth_layers=[23;21;21;19]; %maunally identified 


figure(4)
clf(4)
% set(gcf, 'Position', get(0, 'Screensize'))
set(gcf,'color','w');

t=tiledlayout(2,2);
for ilat = 1:length(lat)
file_nameNPPcafe=['NPP_extracted_CAFElat',num2str(lat(ilat))];
    load([directory_h,file_nameNPPcafe]);
idx(ilat) = calcGlobalWatercolumn(lat(ilat),lon(ilat),sim);
nexttile                              % without this
        plot(sim.ProdNetAnnual(1,idx(ilat).x, idx(ilat).y)*ones(1,12), 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
      
        plot(NPP_extracted(1,:), 'ob--')

ylim([0,3000])
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
mTitle = append('Lat: ', string(lat(ilat)));
title(mTitle)

end
lg  = legend('NUM','CAFE','NumColumns',3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout