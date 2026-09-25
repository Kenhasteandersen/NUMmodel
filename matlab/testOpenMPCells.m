%
% Check and time simulateEulerCells (OpenMP threads over grid cells) against
% the serial loop over cells that simulateGlobal uses without a parallel pool.
%
% Set the number of threads with OMP_NUM_THREADS before starting matlab; the
% number actually used is reported by the timing (there is no way to query it
% from here).
%
% In:
%  nCells - number of cells to test (default 4000)
%  nRep   - number of repetitions for the timing (default 20)
%
% Out:
%  bSuccess - whether the threaded call reproduces the serial loop exactly
%
function bSuccess = testOpenMPCells(nCells, nRep)

arguments
    nCells double = 4000;
    nRep double = 20;
end

lib = loadNUMmodelLibrary();
tEnd = 0.5;
dt = 0.1;

setups = {'setupGeneralistsOnly', 'setupGeneralistsPOM', 'setupNUMmodel'};
bSuccess = true;

for iSetup = 1:length(setups)
    switch setups{iSetup}
        case 'setupGeneralistsOnly', p = setupGeneralistsOnly(10);
        case 'setupGeneralistsPOM',  p = setupGeneralistsPOM(10,1);
        case 'setupNUMmodel',        p = setupNUMmodel();
    end
    %
    % Cells with a range of light, temperature and nutrients:
    %
    k = (1:nCells)';
    L = 300*mod(0.618034*k, 1);
    T = -2 + 32*mod(0.414214*k, 1);
    u0 = (0.5 + mod(0.732051*k, 1))*ones(1,p.n);
    u0(:,p.idxN) = 1 + 150*mod(0.236068*k, 1);
    u0(:,p.idxDOC) = 0.01 + 10*mod(0.3166*k, 1);
    if isfield(p,'idxSi')
        u0(:,p.idxSi) = 10;
    end
    %
    % Reference: serial loop over cells, as in simulateGlobal without a pool
    %
    uRef = u0;
    t0 = tic;
    for iRep = 1:nRep
        for i = 1:nCells
            uRef(i,:) = calllib(lib, 'f_simulateeuler', ...
                uRef(i,:), L(i), T(i), tEnd, dt);
        end
    end
    tRef = toc(t0);
    %
    % Threaded: all cells in one call
    %
    uCells = u0;
    t0 = tic;
    for iRep = 1:nRep
        uCells = calllib(lib, 'f_simulateeulercells', ...
            int32(nCells), uCells, L, T, tEnd, dt);
    end
    tCells = toc(t0);

    dMax = max(abs(uCells(:)-uRef(:)));
    fprintf('%-22s serial %7.3f s   threaded %7.3f s   speedup %5.2fx   max diff %8.1e\n', ...
        setups{iSetup}, tRef, tCells, tRef/tCells, dMax);
    if dMax ~= 0
        bSuccess = false;
        fprintf(2,'  Threaded result differs from the serial loop!\n');
    end
end
