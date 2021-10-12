addpath '..'

nBins = 25; % Use 10 bins for a faster simulation

% The first parameters is the number of size bins. The last parameter is
% the upper size:
p = parametersGeneralistsOnly( nBins, 10.);
pp = p;
p = parametersChemostat( p );

p.mortHTLm = 0*p.mortHTLm; % No HTL mortality

d = [0.01 0.1]; % The mixing rates to run over
%
% Set up figure:
%
clf
tiledlayout(3,3)
%
% Case one: phagotrophy included
%
sweep(p,d, false)
panelSpectra(p, false);
%
% Case two: no phagotrophy
%
p.AF = 0*p.AF; % Setting the affinity for feeding to zero
sweep(p,d, false)
panelSpectra(p, false);
%
% Case three: no phagotrophy or phototrophy for large cells
%
p.pGeneralists.ALm( p.m(3:end)>1e-6 ) = 0; % No phototrophy for larger cells
p.AF = pp.AF;
p.AF( p.m>1e-6 ) = 0; % No phagotrophy for larger cells
sweep(p,d, true)

panelSpectra(p, true);

%
% Sweep over deep nutrient concentrations:
%
function sweep(pp, d, bLastrow)

N0 = logspace(log10(0.01),log10(20),10); % nutrient conditions to sweep over

for j = 1:length(d)
    p = pp;
    p.d = d(j);
    nexttile
    B = zeros(length(N0), length(p.m)-2);
    Light_vs_feeding = B;
    
    p.u0(3:end) = 0.0001; % for faster convergence
    for i = 1:length(N0)
        p.u0(1) = N0(i);
        
        %if (N0(i) < 0.1)
        %    p.tEnd = 5000; % Need to run long for low concentrations
        %else
            p.tEnd = 1000;
        %end
        
        sim(i) = simulateChemostat( p );
        
        B(i,:) = mean(sim(i).B(floor(3*length(sim(i).t)/4):end,:),1);
        tmp = sim(i).rates.JLreal./sim(i).rates.JFreal;
        Light_vs_feeding(i,:) = tmp(3:end);
        
        %plotChemostat(sim)
        %drawnow
    end
    B(B<0) = 0;
    
    Light_vs_feeding(Light_vs_feeding<1) = 1;
    Light_vs_feeding(Light_vs_feeding~=1) = 0;
    
    %%
    % Make the sweep plots:
    %
    surface(p.m(3:end), N0, log10(B));
    shading flat
    axis tight
    set(gca,'xscale','log','yscale','log','XTick',10.^(-8:2))
    
    if bLastrow
        cbar = colorbar('SouthOutside');
        cbar.Limits=[-4,2];
        xlabel('Cell mass ({\mu}g_C)')
    else
        set(gca,'xticklabel','')
    end
    
    %title('log_{10} biomass');
    caxis([-5, 3])
    
    hold on
    %contour3(p.m(3:end), N0, 1000*Light_vs_feeding, [500 500],'w','linewidth',3)
    
    if (j==1)
        ylabel('Deep nutrient concentration ({\mu}g_P/l)')
    end
    
    
end



end
%
% Make a plot of three spectra:
%
function panelSpectra(p, bLastrow)
nexttile

N0 = [0.0669 0.4472 2]; %[0.2 2 20];

for i = 1:length(N0)
    p.u0(1) = N0(i);
    sim = simulateChemostat(p);
    B(i,:) = sim.B(end,:);
    loglog(p.m(3:end), B(i,:), 'b-','linewidth',i);
    hold on
end
set(gca,'xscale','log','yscale','log','XTick',10.^(-8:2))
ylim([1 1000])
xlim([min(p.m) max(p.m)])

if bLastrow
    xlabel('Cell mass ({\mu}g_C)')
else
    set(gca,'xticklabel','')
end


end