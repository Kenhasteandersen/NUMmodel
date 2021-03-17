%
% Make a basic run of the chemostat model
% In:
%  mAdult is the adult sizes of copepods (can be left empty to simulate only
%         unicellular organisms) (default = [], ie only generalists).
%  bUseFortran: Flag indicating whether the fortran library is used
%  (default to false)
%
% Out:
%  sim: Structure holding the results of the simulation
%
function sim = baserunChemostat(mAdult, bUseFortran)

arguments
    mAdult double = []
    bUseFortran logical = false
end
    
if (nargin==0)
    mAdult = [];
end
if (nargin < 2)
    bUseFortran = false;
end
%
% Set parameters:
%
p = parameters(mAdult);
p = parametersChemostat(p);
p.tEnd = 365;
p.bUseLibrary = bUseFortran;
%s
% Setup fortran library:
%
if bUseFortran
    loadNUMmodelLibrary();
    calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult);
end
%
% Simulate
%
tic
sim = simulateChemostat(p, 100);
toc
%
% Plot
%
plotChemostat(sim)