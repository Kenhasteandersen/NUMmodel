%
% load a simulation to be used for initial conditions
%
function sim = loadGlobal(p)

arguments
    p struct
end

load( strcat(fileparts(mfilename('fullpath')),'/',p.pathInit),'sim');
