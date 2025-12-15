load basepath

gridFile=fullfile(base_path,'grid');

load(gridFile,'nx','ny','x','y','bathy')
[X,Y]=ndgrid(x,y);

bathy(bathy==0)=NaN;
bathy(~isnan(bathy))=1;
bathys=bathy(:,:,1);

ko=find(~isnan(bathys));
Xb=X(ko);
Yb=Y(ko);

% fice [fraction]
fice=read_binary('fice.bin',[nx ny 12],'float32');

% wind [m/s]
wind=read_binary('tren_speed.bin',[nx ny 12],'float32');

% dust
fld=double(ncread('iron_1m_ORCA1_MAHO.nc','dust'));
lon=double(ncread('iron_1m_ORCA1_MAHO.nc','nav_lon'));
lat=double(ncread('iron_1m_ORCA1_MAHO.nc','nav_lat'));
mask=double(ncread('ccd_ocal_nemo.nc','OCAL_CCD'));
mask(~isnan(mask))=1;
nt=size(fld,3);
for it=1:nt
  fld(:,:,it)=fld(:,:,it).*mask;
end
ii=find(lon<0);   
lon(ii)=lon(ii)+360;
fldb=interp_2dfield(fld,lon,lat,Xb,Yb);
nt=size(fldb,2);
dust=zeros([nx ny nt]);
fldtmp=zeros([nx ny]);
for it=1:nt
  fldtmp(ko)=fldb(:,it);
  dust(:,:,it)=fldtmp;
end

% qsr [W/m^2]
fn='/data2/spk/OceanCarbon/CORE/ncar_rad.15JUNE2009.nc';
fldb=load_core_variable(fn,'SWDN_MOD',Xb,Yb,1);
nt=size(fldb,2);
qsr=zeros([nx ny nt]);
fldtmp=zeros([nx ny]);
for it=1:nt
  fldtmp(ko)=fldb(:,it);
  qsr(:,:,it)=fldtmp;
end

% hmld [m]
fld=double(ncread('mld_DR003_c1m_reg2.0.nc','mld'));
lon=double(ncread('mld_DR003_c1m_reg2.0.nc','lon'));
lat=double(ncread('mld_DR003_c1m_reg2.0.nc','lat'));
mask=double(ncread('mld_DR003_c1m_reg2.0.nc','mask'));
mask(mask==0)=NaN;
nt=size(fld,3);
for it=1:nt
  fld(:,:,it)=fld(:,:,it).*mask;
end
fld(fld<0)=NaN; % missing data  
fldb=interp_2dfield(fld,lon,lat,Xb,Yb);
nt=size(fldb,2);
hmld=zeros([nx ny nt]);
fldtmp=zeros([nx ny]);
for it=1:nt
  fldtmp(ko)=fldb(:,it);
  hmld(:,:,it)=fldtmp;
end

% ocal_ccd [m]
fld=double(ncread('ccd_ocal_nemo.nc','OCAL_CCD'));
lon=double(ncread('ccd_ocal_nemo.nc','nav_lon'));
lat=double(ncread('ccd_ocal_nemo.nc','nav_lat'));
ii=find(lon<0);   
lon(ii)=lon(ii)+360;
fldb=interp_2dfield(fld,lon,lat,Xb,Yb);
nt=size(fldb,2);
ocal_ccd=zeros([nx ny nt]);
fldtmp=zeros([nx ny]);
for it=1:nt
  fldtmp(ko)=fldb(:,it);
  ocal_ccd(:,:,it)=fldtmp;
end
ocal_ccd=squeeze(ocal_ccd);

save MEDUSA_input_data fice wind dust qsr hmld ocal_ccd
