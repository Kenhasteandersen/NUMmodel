%
% Make a set of basic plots of a simulation
%
function plotSimulation(sim, options)

arguments
    sim = struct;
    % Options for global plots:
    options.sProjection = 'fast';
    options.lat = 60;  % Latitude for watercolumn plots in global simulations
    options.lon = -15;  % Longitude for watercolumn plots in global simulations
    options.tDayPlot = -170; % Days before last time to make size spectrum plot etc.
end

if (options.tDayPlot < 0)
    options.tDayPlot = sim.t(end) + options.tDayPlot;
end

switch sim.p.nameModel
    
    case 'chemostat'
        figure
        clf
        plotGroupsTime(sim);

        figure
        clf
        plotSizespectrum(sim, options.tDayPlot);

        if ~any(isnan(sim.p.seasonalOptions.lat_lon)) || sim.p.seasonalOptions.seasonalAmplitude ~= 0
            figure(3)
            plotSizespectrumTime(sim)
        end
        
    case 'watercolumn'
                
        figure(1)
        clf
        plotWatercolumnTime(sim,'depthMax',200);
        
        figure(2)
        clf
        plotWatercolumn(sim,options.tDayPlot,'depthMax',200);
        
        figure(3)
        % Find the depth of maximum biomass:
        Bdepth = sum( sim.B(options.tDayPlot,:,:),3 );
        iDepth = find(Bdepth==max(Bdepth));
         
        plotSizespectrum(sim,options.tDayPlot,iDepth);
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
        plotGlobal(sim,0,sProjection=options.sProjection);
                
        figure(2)
        clf
        plotWatercolumnTime(sim,options.lat,options.lon, depthMax=200);
        sgtitle( sprintf('Water column at %i, %i',[options.lat,options.lon]))
        
        figure(3)
        tDay = options.tDayPlot;
        if (tDay<1)
            tDay = max(sim.t);
        end
        plotWatercolumn(sim, tDay,options.lat,options.lon, bNewplot=true, depthMax=200);
        sgtitle( sprintf('Watercolumn at (%i,%i) year %i, day %i',[options.lat,options.lon,floor(tDay/365)+1, mod(tDay,365)]));

        figure(4)
        plotSizespectrumTime(sim,1,options.lat,options.lon);
        sgtitle(sprintf('Size spectrum at (%3.0f,%3.0f).\n',[options.lat,options.lon]));
        
        figure(5)
        plotSizespectrum(sim,tDay,1,options.lat,options.lon);
        sgtitle( sprintf('Surface spectrum at (%i,%i) year %i, day %i',[options.lat,options.lon,floor(tDay/365)+1, mod(tDay,365)]));

        figure(6)
        plotGlobalTransect(sim,bAMTtrack=true);
 
    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);
end
