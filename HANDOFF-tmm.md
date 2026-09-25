# Coupling NUM to the TMM (PETSc) driver — handoff

Branch `tmm-coupling`, based on `Develop`. This is a **spike**: it proves the coupling
builds, links and starts. It cannot yet run a simulation — see "Status" below. Nothing
here changes the matlab path, which is untouched and still the way to run the model.

## What this is

`simulateGlobal` does the transport with sparse matrix multiplies in matlab and calls
the fortran library once per grid cell. This branch instead hands NUM to the
[TMM driver](https://github.com/samarkhatiwala/tmm) (Khatiwala), which is C/PETSc and
does transport, time stepping and I/O itself, with MPI.

Worth being clear on one point: **you are already using Khatiwala's transport
matrices** — `parametersGlobal.m` points at his `kelvin.earth.ox.ac.uk` configs. This
does not change the matrices. It replaces the driver around them.

## Why bother

- **Newton-Krylov spinup.** Solving for the periodic steady state instead of
  integrating millennia. For a model that needs spun-up global states this is a
  different order of problem, and it is the main reason to do this at all.
- **The serial floor.** The OpenMP measurements (`HANDOFF-openmp.md`) put a ~2.85 s
  serial floor on a 10-day global run — that is matlab doing the sparse transport
  multiply, and threading cannot touch it. PETSc parallelises transport across ranks,
  so this is the thing that actually removes it.
- **MPI across nodes**, rather than threads on one box.
- **Incremental I/O**, which retires the memory ceiling in `HANDOFF-openmp.md`
  (`simulateGlobal` accumulates the whole time series in RAM: ~10 GB for a year of
  daily output).
- No matlab licence needed to *run*.

## Status

Works:

- `tmm/external_forcing_num.c` implements all five TMM callbacks
  (`iniExternalForcing`, `calcExternalForcing`, `writeExternalForcing`,
  `finalizeExternalForcing`, `reInitializeExternalForcing`).
- It compiles with no warnings, links against `libtmm`, PETSc and
  `libNUMmodel_matlab`, and the resulting `numtmm` binary initialises PETSc and TMM.
- Running it without input gets to TMM's own
  `Must indicate maximum number of steps with the -max_steps option`, i.e. it is live
  and wants runtime input, not fixing.

Does not work yet — in the order you would tackle it:

1. **Input generation.** The transport matrices have to be converted to PETSc binary
   format, plus `dz.petsc` and `latitude.bin`. The `.m` scripts in the TMM repo
   (`models/*/input/make_input_files_*.m`) call `writePetscBin`, which is *not* in that
   repo — it is in Khatiwala's separate matlab utilities. The supported replacement is
   `pytmmutils` (`pip install pytmmutils`), which does the same job in python. This is
   the next real chunk of work and it gates everything else.
2. **Sinking.** Not in the fortran at all. `simulateGlobal` builds sparse upwind
   `Asink` operators in matlab from `f_getsinking` velocities and applies them each
   transport step. A full port needs that as a PETSc operator. The spike sidesteps it
   by using `setupGeneralistsOnly`, which has no sinking species — do not assume it
   generalises to `setupNUMmodel` or anything with POM.
3. **Temperature** is a constant stub (`Tconst`, 10 C). It needs wiring to TMM's
   `Theta` forcing, the way MOPS reads `localTs`.
4. **Analysis.** The largest cost and not a coupling problem: `plotGlobal*`,
   `calcFunction`, `checkConservation` are all built on the matlab `sim` struct. Under
   TMM they would read TMM output instead. Decide whether that is worth it *before*
   investing in 1–3.

## Source data — you already have it

TMM's input generators read a `TransportMatrixConfigs/MITgcm_2.8deg` tree. All of it is
present under `TMs/MITgcm_2.8deg`, verified:

| needed | have |
|---|---|
| `config_data.mat` | yes |
| `grid.mat` | yes |
| `Matrix5/Data/boxes.mat` | yes |
| `Matrix5/Data/profile_data.mat` | yes |
| `GCM/` | yes |

Nothing needs re-downloading.

## Build

PETSc and TMM live in a venv, kept out of the system. MPI comes from Homebrew — **not**
from pip: pip's MPI ships no fortran compiler wrappers, and NUM is fortran.

```sh
brew install mpich

python3 -m venv ~/Documents/tmm-venv
source ~/Documents/tmm-venv/bin/activate
pip install --upgrade pip setuptools wheel cython numpy

export MPICC=$(command -v mpicc)
export MPIF90=$(command -v mpif90)
export MPIFC=$MPIF90
export PETSC_CONFIGURE_OPTIONS='--with-cxx=0 --COPTFLAGS=-O2 --FOPTFLAGS=-O2 --with-debugging=0'
pip install tmmlib pytmmutils
```

Then build NUM's library and the coupled driver:

```sh
cd NUMmodel
cmake -S . -B build && cmake --build build -j
git clone https://github.com/samarkhatiwala/tmm.git ~/Documents/tmm   # for insolation.F
make -C tmm TMM_SRC_REPO=~/Documents/tmm
```

Verified working combination: MPICH 4.x, PETSc 3.25.5, tmmlib 3.2.1, gfortran 16.2.0
(Homebrew), macOS on M2 Ultra.

## Traps, each of which cost a cycle

- **pip's PETSc rejects `--with-cc`/`--with-fc`** in `PETSC_CONFIGURE_OPTIONS`. It
  fails with `RuntimeError: Do not use --with-cc, use the environmental variable
  MPICC`. Set `MPICC`/`MPIF90` in the environment instead, as above.
- **pip's MPI has no fortran wrappers.** Install MPI from brew or source, or fortran
  models cannot be built at all. This is called out in the TMM README and is easy to
  miss.
- **macOS may pop a security dialog** during the PETSc MPI check. If it is not clicked,
  PETSc installs with MPI *silently disabled*. Verify afterwards:
  ```sh
  grep -E "PETSC_HAVE_MPICH|MPIUNI" $PETSC_DIR/include/petscconf.h
  ```
  You want `PETSC_HAVE_MPICH 1` and no `MPIUNI`.
- **The `f_*` symbols are in `libNUMmodel_matlab`, not `libNUMmodel`.** The wrapper
  `NUMmodel_wrap_colmajor.f90` is only in the former. Linking `-lNUMmodel` gives
  undefined `_f_calcderivatives` and friends.
- **`state->c` is an array of `Vec`**, one per tracer, not a single `Vec`. Use
  `state->c[0]` where a template vector is wanted.
- **`localdz` is not a TMM global.** Models declare it themselves and load it with
  `VecLoadVecIntoArray(..., "dz.petsc", ...)`, as MOPS does.
- **`tmm/Makefile` is force-added to git.** `.gitignore` has a blanket `Makefile` rule
  for the CMake in-source build output, which also swallows this one.

## Design note: tendencies, not integration

NUM exposes both `f_simulateeuler` (integrates a cell forward) and
`f_calcderivatives` (returns `dudt`). The coupling uses **`f_calcderivatives`**, so
PETSc owns the time stepping.

This is deliberate. Integrating inside the callback would work and would look more like
the matlab path, but it hides the model's time stepping from PETSc and gives up the
Newton-Krylov spinup — which is the main reason to be doing this. MOPS integrates
internally and converts to a tendency with `mops_biogeochem_copy_data_`; we do not need
that indirection.

One consequence: `f_calcderivatives` is called per grid cell, and NUM's per-cell
calculation is 0-D. TMM hands you profiles (`lNumProfiles`, `lProfileLength`,
`lStartIndices`), so the loop walks columns and then layers. That structure is already
what you want for item 2 (sinking), which is inherently a column operation.
