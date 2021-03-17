%
% Make a basic run of the global transport-matrix model.
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_2.8deg), and be put into the location '../TMs'
%
% In:
%  With no arguments it runs the simple generalist model
%  With a parameter argument in runs the model specified in the parameters.
% Out:
%  A simulation structure
%
function sim = baserunGlobal(p)
%
% Setup a basic run of the global model
%
if (nargin==0)
    p = parameters([]);
    p = parametersGlobal(p); % Use standard low-res model
    %p = parametersGlobal(10,2); % Use MITgcm_ECCO
    p.tEnd = 365;
end

if exist(strcat(p.pathInit,'.mat'), 'file')
    % Load decent initial conditions
    disp('Loading initial conditions from file');
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