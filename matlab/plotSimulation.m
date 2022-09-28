%
% Make a set of basic plots of a simulation
%
function plotSimulation(sim)

switch sim.p.nameModel
    
    case 'chemostat'
        figure(1)
        clf
        plotSizespectrum(sim);
        
        figure(2)
        plotGroupsTime(sim);

        if ~isnan(sim.p.seasonalOptions.lat_lon) | sim.p.seasonalOptions.seasonalAmplitude ~= 0
            figure(3)
            plotSizespectrumTime(sim)
        end
        
    case 'watercolumn'
        day = sim.p.tEnd - 170;
        iDepth = 4;
        
        figure(1)
        clf
        plotWatercolumnTime(sim,'depthMax',200);
        
        figure(2)
        clf
        plotWatercolumn(sim,day,'depthMax',200);
        
        figure(3)
        % Find the depth of maximum biomass:
        Bdepth = sum(sum(sim.B,3),2);
        iDepth = find(Bdepth==max(Bdepth));

        plotSizespectrum(sim,day,iDepth);     

        figure(4)
        plotSizespectrumTime(sim,iDepth);
         
        figure(5)
        plotWatercolumnCommunity(sim)

        figure(6)
        plotWatercolumnCommunity(sim, day)
   
    case 'global'
        figure(1)
        clf
        plotGlobal(sim);
        
        lat = 60;
        lon = -15;
        
        figure(2)
        clf
        plotWatercolumnTime(sim,lat,lon, depthMax=200);
        
        figure(3)
        plotWatercolumn(sim,150,lat,lon, bNewplot=true, depthMax=200);
        
        figure(4)
        plotSizespectrumTime(sim,1,lat,lon);
        title(sprintf('Size spectrum at (%3.0f,%3.0f).\n',[lat,lon]));
        
        figure(5)
        plotSizespectrum(sim,150,1,60,-15);
        
    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);
end
