addpath '..'

nBins = 25; % Use 10 bins for a faster simulation

% The first parameters is the number of size bins. The last parameter is
% the upper size:
p = parametersGeneralistsOnly( nBins, 10.);
p = parametersChemostat( p );

p.bLosses = true; % Allow mixing losses

p.mortHTLm = 0*p.mortHTLm; % No HTL mortality

pp = p;

d = [0.01 0.1]; % The mixing rates to run over
%
% Set up figure:
%
clf
tiledlayout(3,3)
%
% Case one: no phagotrophy
%
p.AF = 0*p.AF; % Setting the affinity for feeding to zero
sweep(p,d, false)
drawnow
panelSpectra(p, d, false);
drawnow
%
% Case two: no phototrophy for large cells
%
p = pp;
p.pGeneralists.ALm( p.m(3:end)>1e-6 ) = 0; % No phototrophy for larger cells
sweep(p,d, false)
drawnow
panelSpectra(p, d, false);
drawnow
%
% Case one: phagotrophy included
%
p = pp;
sweep(p,d, true)
drawnow
panelSpectra(p, d, true);
drawnow


%
% Sweep over deep nutrient concentrations:
%
function sweep(pp, d, bLastrow)

L = 30; % Light level to use

ug_to_uM = 1/30.97; % conversion from ugP/l to Molar

N0 = logspace(log10(0.01),log10(40),10); % nutrient conditions to sweep over

for j = 1:length(d)
    p = pp;
    p.d = d(j);
    p.tEnd = 1000;
    nexttile
    B = zeros(length(N0), length(p.m)-2);
    Light_vs_feeding = B;
    
    p.u0(3:end) = 0.0001; % for faster convergence
    parfor i = 1:length(N0)
        ptmp = p;
        ptmp.u0(1) = N0(i);
        
        %if (N0(i) < 0.1)
        %    p.tEnd = 5000; % Need to run long for low concentrations
        %else
        
        %end
        
        sim = simulateChemostat( ptmp, L );
        
        B(i,:) = mean(sim.B(floor(3*length(sim.t)/4):end,:),1);
        tmp = sim.rates.JLreal./sim.rates.JFreal;
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
    B( B<1e-5 ) = 1e-5;
    surface(p.m(3:end), N0 * ug_to_uM, log10(B));
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
        ylabel('Deep nutrient concentration ({\mu}M)')
    end
    
    
end



end
%
% Make a plot of spectra at three nutrient levels and two mixing rates:
%
function panelSpectra(p, d, bLastrow)
nexttile
uM_to_ug = 30.97; % conversion from Molar to ugP/l.
P0 = [0.006 0.06 0.6] * uM_to_ug; %Three selected deep nutrient levels
%N0 = [0.0669 0.4472 2]; 

sLines = {'b-','r-'}; % Colour of the lines for low and high mixing
for j = 1:length(d)
    p.d = d(j);
    for i = 1:length(P0)
        p.u0(1) = P0(i);
        sim = simulateChemostat(p);
        B(i,:) = mean(sim.B(floor(length(sim.t)/2):end,:),1);
        semilogx(p.m(3:end), B(i,:), sLines{j},'linewidth',i);
        hold on
    end
end

set(gca,'xscale','log','XTick',10.^(-8:2))
ylim([0 100])
xlim([min(p.m) max(p.m)])
ylabel('Biomass ({\mu}g_C/l)')

if bLastrow
    xlabel('Cell mass ({\mu}g_C)')
else
    set(gca,'xticklabel','')
end


end