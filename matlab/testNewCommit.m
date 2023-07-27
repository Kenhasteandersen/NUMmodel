%
% Compare the simulation from the new commit to the one from the previous
% commit. After the function is runned, the new simulation is saved over
% the previous simulation.
%
% The new simulation is runned within the function in chemostat,
% watercolumn or global with the set up setupNUMmodel.
%
% The results of the new simulation are saved in an excel file with the
% comment in input.
%
% The entire new simulation structure can be saved if the path is specified
% in the input saveSimulation, the name of the .mat file will be the date when 
% the function is runned.
%
% The new simulation can be compared with an other simulation if the entire
% path is specified in loadSimulation (path+file name). But the parameters
% of both simulations must be exactly the same. The simulation in the loaded .mat
% file must be named 'simOld'. It is possible to load simulation saved
% thanks to the option saveSimulation.
%
% If the options saveSimulation or loadSimulation are used, the new simuation is not saved over the previous one. 
%
% In:
%   nameModel - name of the model to run ('chemostat','watercolumn' or
%   'global')
%   comment - string, the comment to save in the excel file with the
%   simulation results (empty by default)
%   Optional:
%   options.saveSimulation - path where the new simulation should be saved
%   (the name of the file is the date)
%   option.loadSimulation - path+file name of the simulation to be compared
%   with the new simulation
%
% Out:
%   simNew - new simulation
%   sim - compared simulation, variation between the old and the new
%   simulations in %

function [simNew, sim] = testNewCommit(nameModel,comment,options)

arguments 
    nameModel string = 'global'
    comment string = ''
    options.saveSimulation string = ''; %path
    options.loadSimulation string = ''; % entire path with file's name, the simulation to load must be called simOld
end

%
% Run the new simulation
%
p = setupNUMmodel;
switch nameModel
    case 'chemostat'
        p = parametersChemostat(p);
        p.tEnd = 365;
        p.d = 0.1;
        simNew = simulateChemostat(p);
    case 'watercolumn'
        p = parametersWatercolumn(p);
        p.tEnd = 365;
        simNew = simulateWatercolumn(p,60,-40);
    case 'global'
        p = parametersGlobal(p);
        p.tEnd = 30;
        simNew = simulateGlobal(p);
        simNew.B(simNew.B<0)=0; % Get rid of negative biomasses
end
rates = calcRates(simNew);

%Save new results in an excel file
saveResultsSimulation(simNew,rates,comment);

%
% Compare Simulation
%
if ~isempty(options.loadSimulation)
    if ~exist(options.loadSimulation,'file')
        error('Error: Cannot find the simulation file to load')
    else
        load(options.loadSimulation,'simOld');
    end
else
    load(strcat('testNewCommit\oldSimulation_',nameModel,'.mat'),'simOld'); %load the simulation from the previous commit 
end

% compare the 2 simulations
sim=compareSimulations(simNew,simOld,rates=True);
plotComparedSimulations(sim,simNew);


%
% Save new Simulation
%

if isempty(options.loadSimualtion) && isempty(options.saveSimulation)
    simNew.date=datetime;
    simOld=simNew;
    save(strcat('testNewCommit\oldSimulation_',nameModel,'.mat'),'simOld'); %overwrite on the previous simulation to save space
end

%Saves the entire new simulation structure if the path is specified
if ~isempty(options.saveSimulation)
    if ~exist(options.saveSimulation,'dir')
        mkdir(options.saveSimulation);
    end
    simNew.date=datetime;
    simOld = simNew; %to be able to call this simualtion with options.loadSimulation
    save(strcat(options.saveSimulation,'\',datetime('today')),'simOld');
end
