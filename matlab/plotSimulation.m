%θ
% Make a set of basic plots of a simulation
%
function plotSimulation(sim)

switch sim.p.nameModel
    
    case 'chemostat'
        figure
        clf
        plotGroupsTime(sim);

        figure
        clf
        plotSizespectrum(sim);

        if ~any(isnan(sim.p.seasonalOptions.lat_lon)) || sim.p.seasonalOptions.seasonalAmplitude ~= 0
            figure(3)
            plotSizespectrumTime(sim)
        end
        
    case 'watercolumn'
        day = sim.p.tEnd - 170;
                
        figure(1)
        clf
        plotWatercolumnTime(sim,'depthMax',200);
        
        figure(2)
        clf
        plotWatercolumn(sim,day,'depthMax',200);
        
        figure(3)
        % Find the depth of maximum biomass:
        Bdepth = sum( sim.B(day,:,:),3 );
        iDepth = find(Bdepth==max(Bdepth));
         
        plotSizespectrum(sim,day,iDepth);
        % plotSizespectrum(sim,iDepth);

        figure(4)
        % Find the depth of maximum average biomass:
        Bdepth = sum(sum( sim.B(:,:,:),3 ),1);
        iDepth = find(Bdepth==max(Bdepth));
        plotSizespectrumTime(sim,iDepth);
 
        %figure(5)
        %plotWatercolumnCommunity(sim)

        %figure(6)
        %plotWatercolumnCommunity(sim, day)
   
    case 'global'
        figure(1)
        clf
        plotGlobal(sim);
        
        lat = 60;
        lon = -15;
        Lat = [45.01 41.93 38.63 35.3 31.55 28.61 25.05 21.45 17.74 13.81 10.1 6.26 -1.15 -4.985 -8.765 -12.085 -15.035 -21.6 -24.905 -27.2 -30.01 -32.375 -34.58 -37.01 -39.21 -41.37 -43.79 -46.02];
        Lon = [-13.58 -16.02 -18.53 -20.92 -22.42 -24.36 -26.05 -27.75 -28.93 -28.29 -27.29 -26.38 -24.99 -24.97 -24.95 -24.93 -25.03 -25.05 -25.903 -28.25 -31.15 -33.67 -36.09 -38.83 -41.47 -44.01 -46.95 -49.87];
        
        figure(2)
        clf
        plotWatercolumnTime(sim,lat,lon, depthMax=200);
        
        figure(3)
        plotWatercolumn(sim,150,lat,lon, bNewplot=true, depthMax=200);
        
        figure(4)
        plotSizespectrumTime(sim,1,lat,lon);
        title(sprintf('Size spectrum at (%3.0f,%3.0f).\n',[lat,lon]));
        
        figure(5)
        plotSizespectrum(sim,150,1,lat,lon);
        
        figure(6)
        plotGlobalTransect(sim,Lat,Lon,-1);
        sgtitle('Approximate AMT track - average over 1 year')

    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);
end
