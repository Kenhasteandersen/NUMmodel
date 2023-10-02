% Water column diagnostic plots
% for graphs it calls function: 
% stackedBarplot_percentage(BpnmVec,legend_names,xlabel_names,mTitle)
%
 lat=0;
 lon=-172;
 p=setupNUMmodel;
 p=parametersWatercolumn(p);
 p.tEnd= 5*365;
 sim=simulateWatercolumn(p,lat,lon);
 simR=getSimRates(sim);
%%
day=0;%p.tEnd-170; % set this to average over time
depth_layer=1; % set this to 0 for depth-integrated results
Function_of_WatercolumnDiagnostics_and_Plots(lat,lon,depth_layer,day,p,sim,simR)

% exportgraphics(gcf,[append('biomassNUM','_barplots.png')])       
