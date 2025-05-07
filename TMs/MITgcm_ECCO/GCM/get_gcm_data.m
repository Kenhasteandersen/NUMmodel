load basepath

gridFile=fullfile(base_path,'grid');
load(gridFile,'nx','ny','nz')

Tgcm=read_binary(fullfile(base_path,'MIT','ECCO','theta_ecco_1992_2004_monthly_mean_clim.bin'),[nx ny nz 12],'float64');

Sgcm=read_binary(fullfile(base_path,'MIT','ECCO','salt_ecco_1992_2004_monthly_mean_clim.bin'),[nx ny nz 12],'float64');

%  E-P is in m/s
EmPgcm=read_binary(fullfile(base_path,'MIT','ECCO','eminusp_ecco_1992_2004_monthly_mean_clim.bin'),[nx ny 12],'float64');
Srelaxgcm=[];
saltRelaxTimegcm=[];

save Theta_gcm Tgcm
save Salt_gcm Sgcm
save FreshWaterForcing_gcm EmPgcm Srelaxgcm saltRelaxTimegcm
