lat_to_find = 0;
lon_to_find = -172;
name_png = 'sensitivityNPP_0N172W';
% load('NPP_extracted.mat')
load('NPP_extracted_0_172W.mat')

Sensitivity3x3(lat_to_find,lon_to_find,name_png,NPP_extracted)
%%

lat_to_find = -5;
lon_to_find = 5;
name_png = 'sensitivityNPP_5S5E';
% load('NPP_extracted.mat')
load('NPP_extracted_05S05E.mat')

Sensitivity3x3(lat_to_find,lon_to_find,name_png,NPP_extracted)

%%
lat_to_find = 60;
lon_to_find = -40;
name_png = 'sensitivityNPP_60N40W';
% load('NPP_extracted.mat')
load('NPP_extracted_N60W40.mat')

Sensitivity3x3(lat_to_find,lon_to_find,name_png,NPP_extracted)

