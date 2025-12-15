%
% Make a basic run of the global transport-matrix model.
%
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_2.8deg), and put into the location '../TMs'
%
% In:
%  With no arguments it runs the simple generalist model
%  With a parameter argument it runs the setup specified in the parameters.
%
% Out:
%  A simulation structure
%
function sim = baserunGlobal(p)
%
% Setup a basic run of the global model with only generalists
%
if (nargin==0)
%    mAdult = logspace(log10(0.2), log10(10000), 5);
%    n = 10;
%    nZooplanktons = 6;
%    nPOM = 10;
   p=setupNUMmodelzoo(bParallel=true);
    %p = setupGeneralistsOnly(10, true); % Use 10 size groups and parallel execution
    p = parametersGlobal(p); % Use standard low-res model
end
%
% Simulate
%
if exist(strcat(p.pathInit,'.mat'), 'file')
    % Load decent initial conditions
    disp('Loading initial conditions from file');
    
    sim = loadGlobal(p);
    sim = simulateGlobal(p,sim);
else
    sim = simulateGlobal(p);%,sim); % Simulate
end
sim.B(sim.B<0)=0; % Get rid of negative biomasses
%disp('Calculating functions')
%sim = calcFunction(sim); % Calculate functions
%
% Plots:
%
% disp('Plotting')
plotSimulation(sim, sProjection="eckert4")

checkConservation(sim);
%
% CPU-heavy plots:
%
% animateGlobal(sim);
