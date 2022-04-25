%
% Plot Sheldon biomass spectrum. The biomasses are normalised by the log of
% the ratio between upper and lower masses in each bin
%
function panelSpectrum(sim, ixTime)

    arguments
        sim struct;
        ixTime {mustBeInteger} = length(sim.t); % Defaults to last time step
    end
    p = sim.p;
    
    sLegend = {};
    
    for iGroup = 1:p.nGroups
        ix = p.ixStart(iGroup):p.ixEnd(iGroup);
        m = p.m(ix);
        Delta = p.mUpper(ix)./p.mLower(ix);
        ixB = ix-p.idxB+1;
        %
        % Plot background to show oscillations over the last half of the simulation:
        %
        ixAve = find( sim.t > sim.t(end)/2 );
        if (length(sim.t) >1)
            Blower = min( sim.B(ixAve,ixB) )./ log(Delta);
            Bupper = max( sim.B(ixAve,ixB) )./ log(Delta);
            patch([m, m(end:-1:1)], [Blower, Bupper(end:-1:1)] , ...
                p.colGroup{iGroup},...
                'edgecolor','none', 'facealpha',0.15);
        end
        set(gca,'xscale','log','yscale','log')
        hold on
    end
    %
    % Community spectrum:
    %
%     nPoints = 100000;
%     mc = logspace(log10(p.m(3)), log10(p.m(end)), nPoints);
    
    [mc, Bc] = calcCommunitySpectrum(sim.B, sim);
    legendentries(1) = loglog(mc, Bc, 'linewidth', 4.5,'color',[0.7, 0.7, 0.7]);
    sLegend{1} = 'Community spectrum';
    
    %
    % Group spectra:
    %
    for iGroup = 1:p.nGroups
        ix = p.ixStart(iGroup):p.ixEnd(iGroup);
        m = p.m(ix);
        Delta = p.mUpper(ix)./p.mLower(ix);
        ixB = ix-p.idxB+1;
        %
        % Plot the spectrum:
        %
        legendentries(iGroup+1) = ...
            loglog(m, exp( mean( log(sim.B(ixAve, ixB)./log(Delta)),1)), 'linewidth',2,...
            'color',p.colGroup{iGroup});
        loglog(m, sim.B(ixTime, ixB)./log(Delta), ':','linewidth',1,...
            'color',p.colGroup{iGroup})
    
        sLegend{iGroup+1} = p.nameGroup{iGroup};
    end
    ylim([0.1,500])
    xlim(calcXlim(sim.p))
    hold off
    
    xlabel('Mass ({\mu}gC)')
    ylabel('Sheldon biomass ({\mu}gC/L)')
    
    legend(legendentries, sLegend, 'location','northeast','box','off')