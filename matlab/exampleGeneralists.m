%
%Test case for running the generalist size-based plankton model.
%
close all
clear

%%
% First set up parameters for the plankton:
%
param = parameters([]);
% Then add parameters for the chemostat:
param = parametersChemostat(param);
param.bUseFortran = false; % Make sure not to use the Fortran version of the model
% Run the model to create a "sim" structure:
Light = 100;
sim = simulateChemostat(param, Light);
% Plot the results of the model:
plotChemostat(sim)
%
% The first panel shows the biomass concentration of cells as a fucntion of cell size.
% The second panel shows the mass-specific uptakes in units of per day. In this case we see that there are high uptakes of nutrients (blue line), while little feeding (red).
% The third panel shows the mass-specific losses. HTL er mortality inflicted by higher trophic level organisms.
% The last panel shows nutrient and dissolved organic carbon as a function of time. The black line is the total biomass of plankton.

%%
% Plot of photosynthetic vs phagotrophic production:
%
figure
clf
m = sim.p.m;

NetFixationRate = max(0, (sim.rates.JLreal - sim.p.Jresp))./m; % Subtract respiration
FeedingRate = sim.rates.JFreal./m;
g = sim.rates.Jtot./m;
DOCrate = min(g, sim.rates.JDOC./m);

subplot(2,1,1)
semilogx(m, NetFixationRate, 'g',... 
    m, FeedingRate, 'r', ...
    m, DOCrate, 'm',...
    m, g, 'k', ...
    'linewidth',3)
ylabel('rates (day^{-1})')
ylim([0 1.5])
legend({'Net fixation','Feeding','DOC','Total'})
 
subplot(2,1,2)
semilogx(m, NetFixationRate./g, 'g',...
    m,FeedingRate./g,'r',...
    m,DOCrate./g,'m',...
    'linewidth',3)
xlabel('mass (mugC)')
ylabel('Fraction of growth rate')
ylim([0 1])


%%
% We can get the biomass divided into pico (smaller than 2 mum), nano (2-20 mum), 
% and micro plankton in until of mugC/liter:

[sim.Bpico, sim.Bnano, sim.Bmicro]
%%
% We can adjust the diffusion rate with the parameter "d" in the param structure. 
% The default value is 0.5 per day, which is pretty high:
param.d = 0.05;
sim = simulateChemostat(param, Light);
plotChemostat(sim)
%%
% ... and even lower:
param.d = 0.0001;
sim = simulateChemostat(param, Light);
plotChemostat( sim );

%%
% Run across diffusion rates:
%
d = logspace( -4, 0, 10);
for i = 1:length(d)
    param.d = d(i);
    sim = simulateChemostat(param,Light);
    Btotal(i) = sim.Bgroup(end,:); % Total biomass from last timestep
end

clf
semilogx(d*param.pGeneralists.N0, Btotal, 'linewidth',2)
xlabel('Nutrient flow rate (mugN/L/day')
ylabel('Total biomass (mugC/L)')

    
