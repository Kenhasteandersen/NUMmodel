gridFile='../grid';
load(gridFile,'bathy','nx','ny','nz','x','y','z')

aiceFile=fullfile('MOBI_in','aice.nc');
hiceFile=fullfile('MOBI_in','hice.nc');
hsnoFile=fullfile('MOBI_in','hsno.nc');
swradFile=fullfile('MOBI_in','dnswr.nc');
windFile=fullfile('MOBI_in','ws.nc');
feFile=fullfile('MOBI_in','fe.nc');
feadepFile=fullfile('MOBI_in','O_feflux.nc');
fehydrFile=fullfile('MOBI_in','O_fe_hydr.nc');
sgbathyFile=fullfile('MOBI_in','MITgcm_2.8deg_sg_bathy_PI.nc');

uvicGridFile=fullfile('MOBI_in','UVic_ocean_mask');

aice=squeeze(double(ncread(aiceFile,'AICE')));
aice=aice/100; % convert to fraction

hice=squeeze(double(ncread(hiceFile,'HICE')));

hsno=squeeze(double(ncread(hsnoFile,'HSNO')));

wind=squeeze(double(ncread(windFile,'WS')));

% dissolved Fe concentration when not using prognostic iron model
Fe=squeeze(double(ncread(feFile,'FE')));

% data for prognostic iron cycle
% Atmospheric Fe deposition
fld=squeeze(double(ncread(feadepFile,'O_feflux')));
lon=squeeze(double(ncread(feadepFile,'longitude')));
lat=squeeze(double(ncread(feadepFile,'latitude')));
[X,Y]=ndgrid(x,y);
fldb=interp_2dfield(fld,lon,lat,X(:),Y(:));
nt=size(fld,3);
Fe_adep=zeros([nx ny nt]);
fldtmp=zeros([nx ny]);
for it=1:nt
  fldtmp(:)=fldb(:,it);
  Fe_adep(:,:,it)=fldtmp;
end

% Detrital Fe flux
Fe_detr_flux=zeros([nx ny 12]);

% Hydrothermal Fe input
fld=squeeze(double(ncread(fehydrFile,'O_fe_hydr')));
lon=squeeze(double(ncread(fehydrFile,'longitude')));
lat=squeeze(double(ncread(fehydrFile,'latitude')));
dep=squeeze(double(ncread(fehydrFile,'depth')));
uvicGrid=load(uvicGridFile);
[X,Y,Z]=ndgrid(x,y,z);
kn=find(bathy~=0);
Xb=X(kn);
Yb=Y(kn);
Zb=Z(kn);
[lon1,lat1,dep1]=ndgrid(lon,lat,dep);
ko=find(uvicGrid.bathy~=0);
fldb=griddatan([lon1(ko) lat1(ko) dep1(ko)],fld(ko),[Xb Yb Zb],'linear');
% get rid of NaNs  
k1=find(isnan(fldb));
if ~isempty(k1)
  fldb(k1)=griddatan([lon1(ko) lat1(ko) dep1(ko)],fld(ko),[Xb(k1) Yb(k1) Zb(k1)],'nearest');  
end  
if any(isnan(fldb))
  error('ERROR: interpolation did not fill all target points!')
end
Fe_hydr=zeros([nx ny nz]);
Fe_hydr(kn)=fldb;

sgbathy=squeeze(double(ncread(sgbathyFile,'SGBATHY')));

swrad=squeeze(double(ncread(swradFile,'DNSWR')));

save MOBI_input_data aice hice hsno wind Fe Fe_adep Fe_detr_flux Fe_hydr sgbathy swrad
