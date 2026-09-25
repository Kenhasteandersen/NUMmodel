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

## Why this is worth having — measured on the Mac Studio

The original rationale was speed on many cores: `parfor` copies the 22 MB state array
to and from worker *processes* every transport step, while OpenMP threads share it, so
with 16 performance cores the marshalling should become the limit and OpenMP should
pull ahead. On the 4-performance-core M4 it was a wash (16.4 s vs 16.9 s for a 10-day
global run), which four cores explained.

**That has now been measured on the M2 Ultra (16 performance + 8 efficiency cores) and
OpenMP does not pull ahead.** Full NUM model, 10-day global run:

| path | 10-day global run |
|---|---|
| `parfor`, 24 workers | 7.2 s |
| `bOpenMP`, 16 threads | 7.0 s |

A tie, and the same wash as on the M4 for a different reason. `-O2` and the thread
count were both verified first, so there is no missing optimisation behind that number.

The threading itself works fine. The same run on one thread takes 72.0 s, so OpenMP
delivers 10.05x and the global run is ~96 % threadable. Those two numbers decompose
into a ~2.85 s serial floor per 10 simulated days and ~69 s of parallel work. `parfor`
on 24 workers reaches that floor; OpenMP reaches it on 16. So OpenMP is ~1.5x more
core-efficient — it matches `parfor` while leaving the 8 efficiency cores idle — but
the wall-clock result is a tie. Both give ~4.4 min per simulated year.

**So do not keep this branch for speed over `parfor`. Keep it because it takes the
Parallel Computing Toolbox off the critical path.** `loadNUMmodelLibrary` errors
without the toolbox, and `simulateGlobal` then falls back to a serial loop over cells.
For anyone without that licence the comparison is not 7.0 s against 7.2 s:

| path | 10-day global run |
|---|---|
| serial loop (no toolbox) | 72.0 s |
| `bOpenMP`, 16 threads | 7.2 s |

That 10x, available to any install with a gfortran build and no matlab licence beyond
base, is what this branch actually buys. It also drops parpool startup and the memory
of 24 worker processes, which matters because output memory is already the binding
constraint on long runs (see "Watch the memory on long runs").

The cost is the `threadprivate` discipline in "How the thread safety works". That is a
permanent tax on anyone editing the derivative path, and it buys the 10x above — not
any gain over `parfor`. Weigh it accordingly.

If you want the Mac Studio faster than it is now, threading is not the lever: the run
is already ~96 % parallel and bounded by the ~2.85 s serial floor. The flat kernel on
`openmp-cells` is 2x faster than the threaded OO code on the CPU at every thread count
(see "What is not here").

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

Either works, but only **before the first threaded call**: libgomp reads
`OMP_NUM_THREADS` when the first parallel region runs and latches it for the life of
the process. Setting it later in the same matlab session does nothing.

```matlab
setenv('OMP_NUM_THREADS','16');   % from inside matlab, before the first call
```

```sh
OMP_NUM_THREADS=16 matlab        # or from the shell
```

If it is unset, OpenMP uses all logical cores.

**On the M2 Ultra use 16, not 24.** Measured with `testOpenMPCells(4000, 10)`:

| threads | setupGeneralistsOnly | setupGeneralistsPOM | setupNUMmodel |
|---|---|---|---|
| 1 | 1.26x | 1.19x | 1.01x |
| 2 | 2.14x | 2.05x | 1.93x |
| 4 | 4.09x | 4.11x | 3.59x |
| 8 | 6.84x | 6.49x | 6.57x |
| **16** | **10.90x** | **10.97x** | **14.01x** |
| 24 | 9.91x | 10.16x | 11.68x |

Going 16 → 24 costs ~20 % on the full model: the 8 efficiency cores unbalance the
static schedule, as suspected. Scaling across the 16 performance cores is near-linear.

On the M4 (4 performance + 6 efficiency) unset was fastest, but the efficiency cores
contributed little: 1.85x on 2 threads, 3.2x on 4, 4.4x on 8, 4.6x on 10.

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

On the M2 Ultra, all 24 cores, `-O2`:

```
setupGeneralistsOnly   serial   0.500 s   threaded   0.048 s   speedup 10.46x   max diff  0.0e+00
setupGeneralistsPOM    serial   0.575 s   threaded   0.053 s   speedup 10.88x   max diff  0.0e+00
setupNUMmodel          serial   6.867 s   threaded   0.539 s   speedup 12.74x   max diff  0.0e+00
```

Thread sweep — one matlab process per thread count, because the value is latched at the
first parallel region. A `setenv` loop inside one session silently measures the first
thread count six times:

```sh
cd matlab
for th in 1 2 4 8 16 24; do
  echo "--- ${th} threads"
  OMP_NUM_THREADS=$th matlab -nodisplay -batch "testOpenMPCells(4000, 10)"
done
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

`max diff` must be exactly 0. On the M4: parfor 16.9 s, openmp 16.4 s. On the M2 Ultra:
parfor 7.2 s (24 workers), openmp 7.0 s (16 threads), max diff 0.

To separate threading from everything else, run the same `bOpenMP` case at
`OMP_NUM_THREADS=1` and at 16. On the M2 Ultra that is 72.0 s against 7.2 s, which is
where the ~96 % parallel fraction and the ~2.85 s serial floor above come from.

Finally run `testAll` to confirm nothing else changed.

## Expected run times

Full NUM model on MITgcm_2.8deg (52 749 cells, 54 state variables), per simulated year:

| configuration | per simulated year |
|---|---|
| M4, parfor, 10 workers | 5.3 min |
| M2 Ultra, parfor, 24 workers | ~4.4 min (measured) |
| M2 Ultra, bOpenMP, 16 threads | ~4.4 min (measured) |
| M2 Ultra, serial, no toolbox | ~44 min (measured) |

The earlier estimates here were ~3.5 min for `parfor` and ~2.2 min for `bOpenMP`. The
`bOpenMP` one was simply optimistic: `-O2` and the thread count both check out, and the
run is already ~96 % parallel, so there is no missing factor of two to find. The serial
floor is what stops it — see "Why this is worth having".

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
