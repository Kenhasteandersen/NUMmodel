%
% Set the parameters for a chemostat simulation.
%
% In:
%  p - parameter structure from a call to a setupXX function
%
% Out:
%  p - parameter structure with chemostat fields added.
%
function p = parametersChemostat(p)

p.nameModel = 'chemostat';

p.d = 0.1;  % Mixing rate (1/days)
p.tEnd = 365;  % Time to run in days
p.tSave = 1;

p.depthProductiveLayer = 20; % (meters) Only needed for calculation of functions