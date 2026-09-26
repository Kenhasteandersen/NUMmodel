%
% The number of OpenMP threads to use by default for a global run with
% bOpenMP (see simulateGlobal).
%
% All physical cores is the fastest choice, also on the apple silicon machines
% that mix performance and efficiency cores: on a 4+6 core M-series the six
% efficiency cores still give 31 % over using only the four performance ones.
% The performance core count alone is 'sysctl -n hw.perflevel0.physicalcpu' on
% macos, if you want to experiment.
%
% Out:
%  n - number of threads
%
function n = defaultNumThreads()

n = feature('numcores'); % Physical cores
