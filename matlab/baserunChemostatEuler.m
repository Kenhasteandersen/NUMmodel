%
% Make a basic run of the chemostat model
%
function sim = baserunChemostatEuler(mAdult)
arguments
    mAdult double = []
end
%
% Set parameters:
%
p = setupGeneric(mAdult);
p = parametersChemostat(p);
p.tEnd = 365;
%
% Setup fortran library:
%
if libisloaded('NUMmodel')
    unloadlibrary('NUMmodel')
end
loadNUMmodelLibrary();
calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult);
%
% Simulate
%
tic
sim = simulateChemostatEuler(p, 100);
toc
%
% Plot
%
plotSimulation(sim)