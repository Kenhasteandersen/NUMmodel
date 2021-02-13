%
% Make a basic run of the global transport-matrix model.
% In:
%  With no arguments it runs the simple generalist model
%  With a parameter argument in runs the model specified in the parameters.
% Out:
%  A simulation structureedit s
%
function sim = baserunGlobal(p)
if libisloaded('NUMmodel')
    unloadlibrary('NUMmodel')
end
%
% Setup a basic run of the global model
%
if (nargin==0)
    p = parameters([]);
    p = parametersGlobal(p); % Use standard low-res model
%p = parametersGlobal(10,2); % Use MITgcm_ECCO
    p.tEnd = 365;
end

if exist(p.pathInit, 'file')
    % Load decent initial conditions
    load(p.pathInit);
end
%
% Simulate
%
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