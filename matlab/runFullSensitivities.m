% Runs sensitivity analysis for 5 parameters in 4 water columns
% over noYear number of years
noYears= 4;%25;

% Oligotrophic - HOT station
lat_to_find = 24;
lon_to_find = -158;
name_png = 'sensitivityNPP_HOT';
% load('NPP_extracted.mat')
load('NPP_extracted_HOT.mat')

Sensitivity_paramRemin2(lat_to_find,lon_to_find,name_png,NPP_extracted,noYears)
%%
% Equatorial West Pacific
lat_to_find = 0;
lon_to_find = -172;
name_png = 'sensitivityNPP_0N172W';
% load('NPP_extracted.mat')
load('NPP_extracted_0_172W.mat')

Sensitivity_paramRemin2(lat_to_find,lon_to_find,name_png,NPP_extracted,noYears)
%%
% Eutrophic - West Africa upwelling
lat_to_find = -5;
lon_to_find = 5;
name_png = 'sensitivityNPP_5S5E';
% load('NPP_extracted.mat')
load('NPP_extracted_05S05E.mat')

Sensitivity_paramRemin2(lat_to_find,lon_to_find,name_png,NPP_extracted,noYears)
%%
% Seasonal
lat_to_find = 55;
lon_to_find = -40;
name_png = 'sensitivityNPP_55N40W';
% load('NPP_extracted.mat')
load('NPP_extracted_seasonal.mat')

Sensitivity_paramRemin2(lat_to_find,lon_to_find,name_png,NPP_extracted,noYears)

