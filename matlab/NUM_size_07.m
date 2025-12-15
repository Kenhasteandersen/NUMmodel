function [CW_nongelatinous,CW_gelatinous,CW_median_nongel,CW_median_gel]=NUM_size_07()
load base_Zooscan_project_name.mat
dates = [base_Zooscan.Date];       % extrait toutes les dates (datetime)
idx = month(dates) == 7;           % mois >= juillet
data_juillet_et_plus = base_Zooscan(idx);  % sélectionne les éléments

ln_bvmedian=base_Zooscan(1).d1.X;
bvmedian=exp(ln_bvmedian);

bvmedian = array2table(bvmedian);


namepft=base_Zooscan(1).regroupped.Zoo_groups(194:217,1);

for t = 1:85
NBSS = data_juillet_et_plus(t).regroupped.Ybv_Ellipsoid_BV_spectra(:,194:217); 
sizespectre_bv = NBSS.* data_juillet_et_plus(t).tot.X1; 

sizespectre_bv = array2table(sizespectre_bv);

sizespectre_bv = renamevars(sizespectre_bv, {'sizespectre_bv1','sizespectre_bv2','sizespectre_bv3','sizespectre_bv4','sizespectre_bv5','sizespectre_bv6','sizespectre_bv7','sizespectre_bv8','sizespectre_bv9','sizespectre_bv10','sizespectre_bv11','sizespectre_bv12','sizespectre_bv13','sizespectre_bv14','sizespectre_bv15','sizespectre_bv16','sizespectre_bv17','sizespectre_bv18','sizespectre_bv19','sizespectre_bv20','sizespectre_bv21','sizespectre_bv22','sizespectre_bv23','sizespectre_bv24'},namepft); 
gelatinous(:,t) = sizespectre_bv.gelatinous_filter_feeders +sizespectre_bv.gelatinous_carnivorous;
non_gelatinous(:,t) = sizespectre_bv.all_copepoda + sizespectre_bv.large_crustracean;
end

gelatinous_mean = mean(gelatinous,2);
non_gelatinous_mean = mean(non_gelatinous,2);

WW_gelatinous = gelatinous_mean .* 1025 .*10.^-6;
WW_nongelatinous = non_gelatinous_mean .* 1025 .*10.^-6;

CW_gelatinous = WW_gelatinous .* (0.5 ./100).*10.^3;
CW_nongelatinous = WW_nongelatinous .* (9 ./100).*10.^3;
CW_gelatinous(1:147,:) = 0;
CW_nongelatinous(1:147,:) = 0;

bvmedian = table2array(bvmedian);

WW_median_nongel = bvmedian .* 1025 .*10.^-6;
CW_median_nongel = WW_median_nongel .* (9 ./100).*10.^6;

WW_median_gel = bvmedian .* 1025 .*10.^-6;
CW_median_gel = WW_median_gel .* (0.5 ./100).*10.^6;


%plot(CW_median_nongel, CW_nongelatinous,'r')
%hold on
%plot(CW_median_gel, CW_gelatinous,'b')
%xlim([10.^-6,10^3])
end

