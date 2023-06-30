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
        Lat=[54 15.8 -33];
        Lon=[-21 -41 -8.7];
        
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
        plotGlobalTransect(sim,-1,Lat,Lon);

    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);
end
