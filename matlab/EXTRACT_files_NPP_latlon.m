cd 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\NPP_Satellites';

dataSet_names = {'eppley', 'vgpm', 'cafe'};

NPP_extracted = zeros(length(dataSet_names), 12);

for j = 1:length(dataSet_names)

    myFolder = pwd;
    curr_end = append(string(dataSet_names(j)), '*','.hdf');
    myFiles = dir(fullfile(myFolder, curr_end));

    for i = 1:12

        % (lat, lon)
        % first --> rows: +90 / column: -180
        data = hdfread(myFiles(i).name, 'npp');

        lat_step = 180 / 1080;
        lon_step = 360 / 2160;

        lats = 90:-lat_step:-90;
        lons = -180:lon_step:180;

% !!!!!!!! Put lat & lon that you want to find !!!!!!!!!!!!!!!!!


%--------------------------------------------
%        Initialize coordinates
%--------......................--------
        lat_to_find = 55;
        lon_to_find = -40;
%
%--------------------------------------------
       

        [ d, ix_lat ] = min( abs( lats-lat_to_find ));

        [ d, ix_lon ] = min( abs( lons-lon_to_find ));

        NPP_extracted(j, i) = data(ix_lat, ix_lon);

    end
end

% figure

NPP_extracted(NPP_extracted < 0) = 0;
% plot(NPP_extracted(:,:))

%%

figure
hold on
plot(NPP_extracted(1,:), 'ob--')
plot(NPP_extracted(2,:), 'og--')
plot(NPP_extracted(3,:), 'om--')

legend('Eppley model', 'Standard VGPM', 'CAFE')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
title(mTitle)