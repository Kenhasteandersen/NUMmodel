load basepath

gridFile=fullfile(base_path,'grid');
load(gridFile,'nx','ny','nz')

% T/S computed by propagating surface GCM T/S into interior using TMM
load(fullfile(base_path,'PrescribedSeasonalBC_Matrix5_2','prescribed_bc_equilibrium_fields'),'Tm','Sm')
Tbc=Tm;
Sbc=Sm;

% Surface forcing data
Fice=read_binary('fice.bin',[nx ny 12],'float32');

windspeed=read_binary('tren_speed.bin',[nx ny 12],'float32');

silica=read_binary('sillev1.bin',[nx ny 15],'float32'); % NOTE: 15 slices. This is a mistake in creating the file.
silica=silica(:,:,1:12); % slices 13-15 are all zero.

ironflux=read_binary('mah_flux_smooth.bin',[nx ny 12],'float32');

chl=repmat(0,[nx ny 12]);

save ice_fraction Fice
save wind_speed windspeed
save iron_flux ironflux
save silica silica
save chlorophyll chl
save Theta_bc Tbc
save Salt_bc Sbc
