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

load('simSinking5f.mat') %sinkingPOM_sensitivity folder
load('simSinking10f.mat')
load('simSinking30f.mat')
load('simSinking50f.mat')


lat=[0;-5;24;55];
lon=[172;5;-158;-40];
% lat_to_find=lat(ilat);
noYears=5;
depth_layers=[23;21;21;19]; %maunally identified 

%%
time_dim= size(simFken10.ProdNet,1);
time=time_dim-12+1:time_dim;

figure(4)
clf(4)
set(gcf,'color','w');

t=tiledlayout(2,2);
for ilat = 1:length(lat)
    file_nameNPPcafe=['NPP_extracted_CAFElat',num2str(lat(ilat))];
    load([directory_h,file_nameNPPcafe]);
    idx(ilat) = calcGlobalWatercolumn(lat(ilat),lon(ilat),simFken5);
    nexttile                              % without this
        plot(simFken5.ProdNet(time,idx(ilat).x, idx(ilat).y), 'o-r', 'LineWidth', 2)% mgC/m2/day
        hold on
        plot(simFken10.ProdNet(time,idx(ilat).x, idx(ilat).y), 'o-m', 'LineWidth', 2)% mgC/m2/day
        plot(simFken30.ProdNet(time,idx(ilat).x, idx(ilat).y), 'o-g', 'LineWidth', 2)% mgC/m2/day
        plot(simFken50.ProdNet(time,idx(ilat).x, idx(ilat).y), 'o-c', 'LineWidth', 2)% mgC/m2/day

        plot(NPP_extracted(1,:), 'ok--','LineWidth', 2)

    ylim([0 3000])
    xlabel('Time (month)')
    ylabel('NPP (mgC / m^2 /day)')
    mTitle = append('Lat: ', string(lat(ilat)), ', Lon: ', string(lon(ilat)));
    title(mTitle)
    axis tight
end
lg  = legend('v_{sink}=5m/d','v_{sink}=10m/d','v_{sink}=30m/d','v_{sink}=50m/d','CAFE','NumColumns',3); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout
sgtitle('NPP for varying sinkingPOM','FontSize', 20)
lg.FontSize=14;