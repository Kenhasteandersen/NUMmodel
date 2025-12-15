load basepath

gridFile=fullfile(base_path,'grid');
load(gridFile,'nx','ny','nz')

% T/S computed by propagating surface GCM T/S into interior using TMM
load(fullfile(base_path,'PrescribedSeasonalBC_Matrix1_2','Temperature_bc'),'TRmmg')
Tbc=TRmmg;
load(fullfile(base_path,'PrescribedSeasonalBC_Matrix1_2','Salinity_bc'),'TRmmg')
Sbc=TRmmg;

% Surface forcing data
Fice=read_binary('nasa_icefraction_mth-2d.bin',[nx ny 12],'float32');

windspeed=read_binary('tren_speed_mth-2d.bin',[nx ny 12],'float32');

silica=repmat(7.6838e-3,[nx ny 12]); % from dic_ini_forcing.F

ironflux=read_binary('mah_smooth_mth-2d.bin',[nx ny 12],'float32');

chl=repmat(0,[nx ny 12]);

save ice_fraction Fice
save wind_speed windspeed
save iron_flux ironflux
save silica silica
save chlorophyll chl
save Theta_bc Tbc
save Salt_bc Sbc
