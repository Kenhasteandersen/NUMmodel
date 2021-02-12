%Test case for running the generalist size-based plankton model.
close all
clear

%%
% First set up parameters for the plankton:

param = parameters([]);
% Then add parameters for the chemostat:
param = parametersChemostat(param);
param.bUseFortran = false; % Make sure not to use the Fortran version of the model
% Run the model to create a "sim" structure:
Light = 100;
sim = simulateChemostat(param, Light);
% Plot the results of the model:
plotSimulation(sim)
%
% The first panel shows the biomass concentration of cells as a fucntion of cell size.
% The second panel shows the mass-specific uptakes in units of per day. In this case we see that there are high uptakes of nutrients (blue line), while little feeding (red).
% The third panel shows the mass-specific losses. HTL er mortality inflicted by higher trophic level organisms.
% The last panel shows nutrient and dissolved organic carbon as a function of time. The black line is the total biomass of plankton.

%%
% We can get the biomass divided into pico (smaller than 2 mum), nano (2-20 mum), 
% and micro plankton in until of mugC/liter:

[sim.Bpico, sim.Bnano, sim.Bmicro]
%%
% We can adjust the diffusion rate with the parameter "d" in the param structure. 
% The default value is 0.5 per day, which is pretty high:
param.d = 0.05;
sim = simulateChemostat(param, Light);
plotSimulation(sim)
%%
% ... and even lower:
param.d = 0.0001;
sim = simulateChemostat(param, Light);
plotSimulation( sim );

%%
% Run across diffusion rates:

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

    
