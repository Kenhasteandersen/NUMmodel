% addpath ~/Documents/Source/NUMmodel/matlab

% The first parameters is the number of size bins. The last parameter is
% the upper size:
p = parametersGeneralistsOnly( 25, 10.);
p = parametersChemostat( p );

p.mortHTLm = p.mortHTLm / 4; % Lowering the background mortality (not needed, actually).

d = [0.005 0.5]; % The mixing rates to run over
clf
tiledlayout(length(d),2)
sweep(p,d) % With phagotrophy

p.AF = 0*p.AF; % Setting the affinity for feeding to zero
sweep(p,d) % without phagotrophy

%
% Sweep over deep nutrient concentrations:
%
function sweep(pp, d)

N0 = logspace(-1,3,10);

for j = 1:length(d)
    p = pp;
    p.d = d(j);
    nexttile
    B = zeros(length(N0), length(p.m)-2);
    Light_vs_feeding = B;
    
    p.u0(3:end) = 0.0001; % for faster convergence
    for i = 1:length(N0)
        p.u0(1) = N0(i);
        
        if (N0(i) < 0.1)
            p.tEnd = 5000; % Need to run long for low concentrations
        else
            p.tEnd = 365;
        end
        
        sim = simulateChemostat( p );
        
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
    surface(p.m(3:end), N0, log10(B));
    shading flat
    axis tight
    set(gca,'xscale','log','yscale','log','XTick',10.^(-8:2))
    cbar = colorbar('SouthOutside');
    
    cbar.Limits=[-4,2]
    
    title('log_{10} biomass');
    caxis([-5, 3])
    
    
    hold on
    contour3(p.m(3:end), N0, 1000*Light_vs_feeding, [500 500],'w','linewidth',3)
    
    xlabel('Cell mass ({\mu}gC)')
    ylabel('Deep nutrient concentration ({\mu}gN/l)')
end

end