function sim = baserunGlobal(sim)
if libisloaded('NUMmodel')
    unloadlibrary('NUMmodel')
end
%
% Make a basic run of the global model
%
if (nargin==0)
    p = parameters([]);
    p = parametersGlobal(p); % Use standard low-res model
%p = parametersGlobal(10,2); % Use MITgcm_ECCO
    p.tEnd = 365;
else
    p = sim.p;
end

if exist(p.pathInit, 'file');
    % Load decent initial conditions
    load(p.pathInit);
end

sim = simulateGlobal(p);%,sim); % Simulate
sim.B(sim.B<0)=0; % Get rid of negative biomasses
%disp('Calculating functions')
%sim = calcGlobalFunction(sim); % Calculate functions
%
% Plots:
%
disp('Plotting')
figure(1)
clf
plotGlobal(sim);

figure(2)
clf
plotGlobalWatercolumnTime(60,-10,sim);
%
% CPU-heavy plots:
%
% animateGlobal(sim);