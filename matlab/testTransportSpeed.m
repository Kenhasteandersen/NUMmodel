%
% Measure how much time a global run spends on the transport matrix step, and how
% much faster that step could be if it were moved into the fortran library and
% threaded with OpenMP.
%
% It times, for one month of transport matrices:
%  - the biology call, so the share of the step that transport takes is known
%  - matlab's Aimp*(Aexp*u), as simulateGlobal does it now
%  - the same with the boxes reordered for locality
%  - the threaded fortran kernel in Fortran/spmmbench.f90, at a range of thread
%    counts and in three layouts of the state
%
% The fortran part is skipped if spmmbench has not been built. Build it with
%
%   gfortran -O3 -fopenmp -o spmmbench Fortran/spmmbench.f90
%
% from the root of the repository. -O3 is needed: at -O2 the kernel is twice as
% slow because the inner loop does not vectorize.
%
% Set the number of threads with OMP_NUM_THREADS before starting matlab, as for
% the rest of the OpenMP code.
%
% In:
%  p        - parameter structure from one of the setups, e.g. setupNUMmodel().
%             Only p.n (the number of state variables) is used.
%  nTMmodel - 1 for MITgcm_2.8deg, 2 for MITgcm_ECCO (default 2)
%  options.bFortran - whether to run the fortran benchmark (default true)
%  options.dirData  - where to put the exported matrices. Needs about 25 bytes
%                     per nonzero, which is ~3 GB for ECCO (default tempdir)
%
% Out:
%  res - structure with the timings
%
function res = testTransportSpeed(p, nTMmodel, options)

arguments
    p struct = setupNUMmodel();
    nTMmodel double = 2;
    options.bFortran logical = true;
    options.dirData = fullfile(tempdir, 'spmmbench');
end

pTM = parametersGlobal(p, nTMmodel);
nG = pTM.n;
root = fileparts(fileparts(mfilename('fullpath')));
%
% Load one month of transport matrices and apply the same transformations as
% simulateGlobal:
%
fprintf('Loading %s ...\n', [pTM.pathMatrix '01.mat']);
t = tic;
S = load([pTM.pathMatrix '01.mat'], 'Aexp', 'Aimp');
load(pTM.pathGrid, 'deltaT');
nb = size(S.Aexp,1);
Aexp = speye(nb,nb) + (pTM.dtTransport*24*60*60)*S.Aexp;
Aimp = S.Aimp^(pTM.dtTransport*24*60*60/deltaT);
clear S
fprintf('  %.1f s.  nb = %d boxes, nGrid = %d, state = %.0f MB\n', ...
        toc(t), nb, nG, nb*nG*8/2^20);
fprintf('  nnz(Aexp) = %d (%.1f per row), nnz(Aimp) = %d (%.1f per row)\n', ...
        nnz(Aexp), nnz(Aexp)/nb, nnz(Aimp), nnz(Aimp)/nb);

rng(7);
u = rand(nb, nG);
res = struct('nb',nb,'nG',nG);
%
% The biology, for the ratio. This needs the library to be compiled with OpenMP;
% without it the call still works but runs on one thread.
%
lib = loadNUMmodelLibrary();
nThreads = calllib(lib, 'f_getmaxthreads', int32(0));
calllib(lib, 'f_simulateeulercells', int32(nb), u, ones(nb,1)*60, ones(nb,1)*15, ...
        pTM.dtTransport, pTM.dt);
res.biology = best(@() calllib(lib, 'f_simulateeulercells', int32(nb), u, ...
        ones(nb,1)*60, ones(nb,1)*15, pTM.dtTransport, pTM.dt), 3);
%
% Matlab's transport step, as supplied and reordered:
%
res.matlab = best(@() Aimp*(Aexp*u), 3);

fprintf('\nreordering the boxes ...\n');
t = tic; pp = symrcm(Aexp+Aexp.'+Aimp+Aimp.'); fprintf('  symrcm %.1f s\n', toc(t));
AexpP = Aexp(pp,pp); AimpP = Aimp(pp,pp); uP = u(pp,:);
res.matlabReordered = best(@() AimpP*(AexpP*uP), 3);

fprintf('\n');
fprintf('biology, %d threads              : %7.3f s per transport step\n', ...
        nThreads, res.biology);
fprintf('matlab Aimp*(Aexp*u)             : %7.3f s   (%.0f%% of the step)\n', ...
        res.matlab, 100*res.matlab/(res.matlab+res.biology));
fprintf('  with the boxes reordered       : %7.3f s   (%.2fx)\n', ...
        res.matlabReordered, res.matlab/res.matlabReordered);
%
% The fortran kernel:
%
if ~options.bFortran
    return
end
exe = fullfile(root, 'spmmbench');
if ~isfile(exe)
    fprintf(['\nspmmbench not found. Build it to measure the fortran kernel:\n' ...
             '    cd %s\n    gfortran -O3 -fopenmp -o spmmbench Fortran/spmmbench.f90\n'], root);
    return
end
d = options.dirData;
if ~isfolder(d), mkdir(d); end
fprintf('\nexporting the reordered matrices to %s ...\n', d);
t = tic;
writeCSR(fullfile(d,'Aexp'), AexpP);
writeCSR(fullfile(d,'Aimp'), AimpP);
writeBin(fullfile(d,'u.bin'), uP);
writeBin(fullfile(d,'ref.bin'), AimpP*(AexpP*uP));
fid = fopen(fullfile(d,'dims.txt'),'w'); fprintf(fid,'%d %d\n', nb, nG); fclose(fid);
fprintf('  %.1f s\n\n', toc(t));

[st, out] = system(sprintf('"%s" "%s"', exe, d));
disp(out)
if st ~= 0
    warning('spmmbench exited with status %d', st);
end
fprintf(['Compare the tracer-major column against the matlab time above. The\n' ...
         'reordered matlab time is the fair baseline, since the fortran kernel\n' ...
         'is given the reordered matrices too.\n']);
end

function t = best(f, n)
    t = inf;
    for i = 1:n
        t0 = tic; f(); t = min(t, toc(t0));
    end
end

function writeBin(f, a)
    fid = fopen(f,'w'); fwrite(fid, a(:), 'double'); fclose(fid);
end

%
% Write A in compressed sparse row form. Matlab stores sparse matrices by column,
% so the columns of A' are the rows of A.
%
function writeCSR(stem, A)
    [colOfA, rowOfA, val] = find(A.');
    n = size(A,1);
    rowptr = int32([0; cumsum(accumarray(rowOfA, 1, [n 1]))]);
    fid = fopen([stem '_rowptr.bin'],'w'); fwrite(fid, rowptr,        'int32');  fclose(fid);
    fid = fopen([stem '_col.bin'],'w');    fwrite(fid, int32(colOfA), 'int32');  fclose(fid);
    fid = fopen([stem '_val.bin'],'w');    fwrite(fid, val,           'double'); fclose(fid);
    fid = fopen([stem '_nnz.txt'],'w');    fprintf(fid,'%d\n', numel(val));      fclose(fid);
end
