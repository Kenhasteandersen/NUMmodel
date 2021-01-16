function sim = baserunGlobal()
if libisloaded('NUMmodel')
    unloadlibrary('NUMmodel')
end
%
% Make a basic run of the global model
%
p = parameters([]);
p = parametersGlobal(p); % Use standard low-res model
%p = parametersGlobal(10,2); % Use MITgcm_ECCO

%load(p.pathInit); % Load decent initial conditions
%p.tSave = 1; % Save every 10th day
%p.tEnd = 5;

p.tEnd = 10*365;
sim = simulateGlobal(p);%,sim); % Simulate
sim.B(sim.B<0)=0; % Get rid of negative biomasses
disp('Calculating functions')
%sim = calcGlobalFunction(sim); % Calculate functions
%
% Plots:
%
disp('Plotting')
close all
figure
plotGlobal(sim);

figure
plotGlobalWatercolumnTime(60,-10,sim);
%
% CPU-heavy plots:
%
% animateGlobal(sim);