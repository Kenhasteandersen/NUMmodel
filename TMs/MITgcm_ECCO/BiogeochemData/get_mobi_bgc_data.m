aiceFile=fullfile('MOBI_TMM_ECCOgrid_in','aice.nc');
hiceFile=fullfile('MOBI_TMM_ECCOgrid_in','hice.nc');
hsnoFile=fullfile('MOBI_TMM_ECCOgrid_in','hsno.nc');
swradFile=fullfile('MOBI_TMM_ECCOgrid_in','dnswr.nc');
windFile=fullfile('MOBI_TMM_ECCOgrid_in','ws.nc');
feFile=fullfile('MOBI_TMM_ECCOgrid_in','fe.nc');
sgbathyFile=fullfile('MOBI_TMM_ECCOgrid_in','MIT_gcm_ECCO_sg_bathy_PI.nc');

aice=getnc(aiceFile,'AICE',-1,-1,-1,[4 3 2 1]);
aice=aice/100; % convert to fraction

hice=getnc(hiceFile,'HICE',-1,-1,-1,[4 3 2 1]);

hsno=getnc(hsnoFile,'HSNO',-1,-1,-1,[4 3 2 1]);

wind=getnc(windFile,'WS',-1,-1,-1,[4 3 2 1]);

Fe=getnc(feFile,'FE',-1,-1,-1,[4 3 2 1]);

sgbathy=getnc(sgbathyFile,'SGBATHY',-1,-1,-1,[3 2 1]);

swrad=getnc(swradFile,'DNSWR',-1,-1,-1,[4 3 2 1]);

save MOBI_input_data aice hice hsno wind Fe sgbathy swrad
