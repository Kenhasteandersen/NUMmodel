load basepath

gridFile=fullfile(base_path,'grid');
load(gridFile,'nx','ny','nz')

Tgcm=rdmds(fullfile(base_path,'MIT','SeasonalDST3','Ttave'),NaN);

Sgcm=rdmds(fullfile(base_path,'MIT','SeasonalDST3','Stave'),NaN);

%  E-P is in m/s
EmPgcm=[];
Srelaxgcm=read_binary(fullfile(base_path,'MIT','input','lev_monthly_salt.bin'),[nx ny 12],'float32');
saltRelaxTimegcm=62208000.0;

save Theta_gcm Tgcm
save Salt_gcm Sgcm
save FreshWaterForcing_gcm EmPgcm Srelaxgcm saltRelaxTimegcm
