# OpenMP threading over grid cells — handoff for the Mac Studio

Branch `openmp-cells2`, based on `Develop`. Supersedes `HANDOFF-openmp-cells.md`,
which describes the older `openmp-cells` branch (see "What is not here" below).

## What this does

`simulateGlobal` used to call the fortran library once per grid cell, ~52 700 times
per transport step, parallelised with matlab's `parfor`. This branch adds
`simulateEulerCells`, which integrates all cells in **one** call and distributes them
over OpenMP threads inside the library. The object-oriented fortran structure is
unchanged.

Enable it with

```matlab
sim = simulateGlobal(p, bOpenMP=true);
```

Results are bitwise identical to both the serial loop and the `parfor` path.

## Why it matters on the Mac Studio, and not on a laptop

On a 4-performance-core M4 this is a wash against `parfor` (16.4 s vs 16.9 s for a
10-day global run) — four cores is the ceiling either way. The gain comes when there
are many cores: `parfor` copies the 22 MB state array to and from worker *processes*
every transport step, while OpenMP threads share it. With 16 performance cores that
marshalling becomes the limit, so expect OpenMP to pull ahead. Measure it (below)
rather than assuming.

## Build

Needs gfortran with OpenMP. Apple's clang has no OpenMP, but only fortran is
compiled here, so Homebrew gcc is enough:

```sh
brew install gcc
cd NUMmodel
cmake -S . -B build          # should print: Found OpenMP_Fortran: -fopenmp
cmake --build build -j
cmake --install build        # copies the libraries into lib/
```

Check the flags picked up both OpenMP and optimisation:

```sh
grep '^Fortran_FLAGS =' build/Fortran/CMakeFiles/NUMmodel.dir/flags.make
# expect: -fopenmp -O2 -DNDEBUG -fPIC -fopenmp
```

`-O2` matters a great deal: until recently the flags never reached the compiler and
everything was built at `-O0`, which is 4x slower. If you see no `-O2`, stop and
find out why before measuring anything.

To build without threading: `cmake -S . -B build -DNUM_OPENMP=OFF`. The `!$omp`
directives are then comments and everything runs serially.

## Setting the number of threads

Either works — the environment variable is read when the first parallel region runs,
not when matlab starts:

```matlab
setenv('OMP_NUM_THREADS','16');   % from inside matlab, before the first call
```

```sh
OMP_NUM_THREADS=16 matlab        # or from the shell
```

If it is unset, OpenMP uses all logical cores. On the M4 (4 performance + 6
efficiency) unset was fastest, but the efficiency cores contribute little: measured
speedups were 1.85x on 2 threads, 3.2x on 4, 4.4x on 8, 4.6x on 10. **On the M2 Ultra
(16 performance + 8 efficiency) compare 16 threads against 24** — the efficiency
cores may well slow it down by unbalancing the static schedule.

Do not combine `bOpenMP=true` with a `parpool`; it replaces the `parfor` loop. Use
`setupNUMmodel()` rather than `setupNUMmodel(bParallel=true)`.

## Verify, then measure

`testOpenMPCells(nCells, nRep)` runs the threaded call and the serial loop over the
same cells and checks they agree exactly. Run it first — a non-zero `max diff` means
a thread-safety problem, and nothing else is worth measuring until it is fixed.

```matlab
testOpenMPCells          % 4000 cells, 20 repetitions
```

Reference output on the M4, all cores, `-O2`:

```
setupGeneralistsOnly   serial   0.446 s   threaded   0.094 s   speedup  4.77x   max diff  0.0e+00
setupGeneralistsPOM    serial   0.500 s   threaded   0.112 s   speedup  4.46x   max diff  0.0e+00
setupNUMmodel          serial   5.266 s   threaded   1.528 s   speedup  3.44x   max diff  0.0e+00
```

Thread sweep:

```matlab
for th = [1 2 4 8 16 24]
    setenv('OMP_NUM_THREADS', num2str(th));
    fprintf('--- %d threads\n', th); testOpenMPCells(4000, 10);
end
```

Then the real comparison, `parfor` against OpenMP on a global run:

```matlab
p = setupNUMmodel(bParallel=true); p = parametersGlobal(p); p.tEnd = 10; p.tSave = 5;
tic; sP = simulateGlobal(p); tP = toc;
delete(gcp('nocreate'));
p = setupNUMmodel(); p = parametersGlobal(p); p.tEnd = 10; p.tSave = 5;
setenv('OMP_NUM_THREADS','16');
tic; sO = simulateGlobal(p, bOpenMP=true); tO = toc;
fprintf('parfor %.1f s, openmp %.1f s, max diff %.3e\n', ...
    tP, tO, max(abs(double(sP.B(:))-double(sO.B(:)))));
```

`max diff` must be exactly 0. On the M4: parfor 16.9 s, openmp 16.4 s.

Finally run `testAll` to confirm nothing else changed.

## Expected run times

Full NUM model on MITgcm_2.8deg (52 749 cells, 54 state variables), per simulated year:

| configuration | per simulated year |
|---|---|
| M4, parfor, 10 workers | 5.3 min |
| M2 Ultra, parfor (estimate) | ~3.5 min |
| M2 Ultra, bOpenMP (estimate) | ~2.2 min |

If the Ultra lands far off ~2 min/year, something is wrong — check `-O2` first, then
the thread count.

## Watch the memory on long runs

Output, not computation, is the limit. `simulateGlobal` accumulates the whole time
series in memory before returning: `sim.B` alone is 128x64x15x51 singles per save,
so about **27.6 MB per save** including the other fields.

| run | saves | memory |
|---|---|---|
| 1 year, daily | 365 | 10 GB |
| 10 years, monthly | 120 | 3.3 GB |
| 100 years, monthly | 1200 | 33 GB |

With 64 GB a century run with monthly output puts half your RAM in one array. Either
coarsen `p.tSave` or write incrementally to disk.

## How the thread safety works — read this before changing the library

The `group` objects hold the **rates of the cell currently being integrated**, not
just parameters, so they cannot be shared between threads. Therefore:

- `group`, `upositive`, `F` in `NUMmodel.f90` are `threadprivate`
- `fTemp2`, `fTemp15` and the cached `Told` in `globals.f90` are `threadprivate`
- each non-master thread deep-copies `group` from `groupShared`, a snapshot taken by
  the master at the start of the parallel region (this is also what picks up parameter
  changes from `setHTL` between calls)
- `theta`, `pHTL`, `thetaPOM`, `ixStart`, `ixEnd`, `nGrid` are only read after setup,
  so they stay shared

`copyin` was tried instead of the deep copy and segfaulted in gfortran on the
polymorphic allocatable components, hence the manual copy.

**If you add a module-level variable that `calcDerivatives` writes to, it must be
made `threadprivate` too, or results will silently depend on the thread count.** That
is the one way this design breaks quietly. `testOpenMPCells` is the guard: it compares
threaded against serial bitwise, so run it after any change to the derivative path.

Each thread holds a full copy of `group`. That is small, but it grows with thread
count — worth remembering at 24 threads.

## What is not here, and why

The older `openmp-cells` branch also contains a GPU offload prototype:
`NUMmodel_offload.f90` with a flat `!$omp target` kernel, plus a refactor of
`generalists.f90` and `spectrum.f90` into `elemental` "core" procedures to support it.
None of that is on this branch.

Two reasons. First, **the Mac Studio cannot run it**: Apple GPUs have no FP64 at all
(Metal has no double type) and no fortran compiler can emit OpenMP offload code for
them — the targets are nvptx, amdgcn and spir64. The model is `real(dp)` throughout,
and the mass-balance checks now pass at 1e-12, which FP32 could not sustain.
Second, the flat kernel **duplicates the community assembly** — food, predation
mortality, the gamma corrections, HTL routing — from `calcDerivatives`, so it is a
second copy of the model core that has to be kept in sync by hand, and it covers
generalists only.

Worth knowing: the flat kernel is a genuine **2x faster than the threaded OO code on
the CPU**, at every thread count and both optimisation levels. So the flattening pays
off before any GPU is involved. If that 2x is ever worth the duplication, the honest
way to take it is a full flat rewrite with the OO version retained as a reference
oracle and a bitwise test between the two — not a hand-maintained partial copy.

If you do get access to an NVIDIA card (A100/H100, FP64), the offload work on
`openmp-cells` is the starting point, but note that branch predates several fixes on
`Develop` (the `gammaDOC` sign error it lists as a "pre-existing quirk" is now fixed)
and will need merging.

## Open items

- `lib/` still ships `-O0` Windows and Linux binaries, and they predate the
  `f_simulatechemostateuler` signature change. They need rebuilding on those platforms.
- `NUMmodeltest_cells.f90` on `openmp-cells` is a useful standalone benchmark that
  needs no matlab. It was not brought over because it calls the offload kernel; it is
  easy to port by deleting the flat-kernel section.
- `calcRatesGeneralists` overwrites `JF` with `JFreal`, so the corrector pass uses the
  predictor's realised feeding as available food. Noted in the old handoff, never
  investigated, unrelated to threading.
