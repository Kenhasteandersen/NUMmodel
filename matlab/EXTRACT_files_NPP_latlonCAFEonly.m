% load('C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Model EVALUATION\CAFE\NPP_monthlyAVG.mat')
data=nppAvg;
NPP_extracted = zeros(1,12);

lat=zeros(1080,1);
lat(1)=1/6;
for i=2:1080
    lat(i)=lat(i-1)+1/6;
end
lat=lat-90;
%
lon=zeros(2160,1);
lon(1)=1/6;
for i=2:2160
    lon(i)=lon(i-1)+1/6;
end
    for i = 1:12

        lats=lat;
        lons=lon;

% !!!!!!!! Put lat & lon that you want to find !!!!!!!!!!!!!!!!!


%--------------------------------------------
%        Initialize coordinates
%--------......................--------
        lat_to_find = 24;
        lon_to_find = -158;
%
%--------------------------------------------
       

        [ d, ix_lat ] = min( abs( lats-lat_to_find ));

        [ d, ix_lon ] = min( abs( lons-lon_to_find ));

        NPP_extracted(i) = data(i,ix_lat, ix_lon);

    end
% end

% figure

NPP_extracted(NPP_extracted < 0) = 0;
% plot(NPP_extracted(:,:))
save(['NPP_extracted_CAFElat',num2str(lat_to_find),'.mat'],'NPP_extracted'); 

%%

figure
hold on
plot(NPP_extracted(1,:), 'ob--')

legend('CAFE')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
title(mTitle)