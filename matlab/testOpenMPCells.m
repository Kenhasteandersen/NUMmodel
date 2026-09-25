%
% Test and time the integration of many grid cells in one library call:
%   1) reference: loop over cells calling f_simulateeuler (as simulateGlobal without a pool)
%   2) f_simulateeulercells: OpenMP threads over cells (object-oriented code, all setups)
%   3) f_simulateeulercellsgeneralists: flat kernel that can be offloaded to a GPU
%      (generalists only)
% Set the number of threads with OMP_NUM_THREADS before starting matlab.
% Run without a parallel pool.
%
% In:
%  nCells: number of grid cells (default 10000)
%
% Out:
%  bSuccess: true if both library calls reproduce the reference
%
function bSuccess = testOpenMPCells(nCells)
arguments
    nCells = 10000;
end

p = setupGeneralistsOnly(10);
sLibname = loadNUMmodelLibrary();

k = (1:nCells)';
L = 300 * mod(0.618034*k, 1);
T = -2 + 32 * mod(0.414214*k, 1);
u0 = 0.5 + mod(0.732051*k, 1) * ones(1,p.n);
u0(:,p.idxN) = 1 + 150*mod(0.236068*k, 1);
u0(:,p.idxDOC) = 0.01 + 10*mod(0.3166*k, 1);

tEnd = 0.5; % As dtTransport in simulateGlobal
dt = 0.1;

uRef = u0;
tic
for i = 1:nCells
    uRef(i,:) = calllib(sLibname, 'f_simulateeuler', uRef(i,:), L(i), T(i), tEnd, dt);
end
fprintf('Loop over cells       : %6.3f s\n', toc);

tic
uCells = calllib(sLibname, 'f_simulateeulercells', int32(nCells), u0, L, T, tEnd, dt);
fprintf('f_simulateeulercells  : %6.3f s, max diff %g\n', toc, max(abs(uCells(:)-uRef(:))));

tic
uFlat = calllib(sLibname, 'f_simulateeulercellsgeneralists', int32(nCells), u0, L, T, tEnd, dt);
fprintf('flat (offload) kernel : %6.3f s, max diff %g\n', toc, max(abs(uFlat(:)-uRef(:))));

bSuccess = isequal(size(uCells), size(uRef)) && max(abs(uCells(:)-uRef(:))) < 1e-10 ...
    && max(abs(uFlat(:)-uRef(:))) < 1e-10;
