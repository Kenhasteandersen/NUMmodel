% Water column diagnostic plots
% for graphs it calls function: 
% stackedBarplot_percentage(BpnmVec,legend_names,xlabel_names,mTitle)
%
 lat=60;
 lon=-40;
 p=setupNUMmodel;
 p=parametersWatercolumn(p);
 newSinkingPOM = 0.78622; % optional
 setSinkingPOM(p, newSinkingPOM);
 p.tEnd= 15*365;
 sim=simulateWatercolumn(p,lat,lon);
 simR=getSimRates(sim);

% Calculate Biomass of organism that rely more than 60% on photosynthesis
 Bph=calcPhyto(sim,simR);


 % Consider adding NPP
%%
day=10718;%p.tEnd-170; % set this to average over time
depth_layer=1; % set this to 0 for depth-integrated results
Function_of_WatercolumnDiagnostics_and_Plots(lat,lon,depth_layer,day,p,sim,simR)

exportgraphics(gcf,[append('biomassNUM','_barplots.png')])       
