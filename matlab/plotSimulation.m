%
% Make a set of basic plots of a simulation
%
function plotSimulation(sim)

switch sim.p.nameModel
    
    case 'chemostat'
        clf
        plotSizespectrum(sim);
        
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
        
        figure(2)
        clf
        plotWatercolumnTime(sim,60,-10, depthMax=200);
        
        figure(3)
        plotWatercolumn(sim,150,60,-10, bNewplot=true, depthMax=200);
        
        figure(4)
        plotSizespectrumTime(sim,1,60,-10);
        
        figure(5)
        plotSizespectrum(sim,150,1,60,-10);
        
    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);
end
