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
        
    case 'watercolumn'
        day = 170;
        
        figure(1)
        plotWatercolumnTime(sim,'depthMax',200);
        
        figure(2)
        plotWatercolumn(sim,day,'depthMax',200);
        
        figure(3)
        plotSizespectrum(sim,day,1);
        
        figure(4)
        plotSizespectrumTime(sim,1);
        
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
