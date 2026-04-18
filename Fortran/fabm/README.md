# NUMmodel FABM coupling

This directory contains the [FABM](https://github.com/fabm-model/fabm) institute
library that couples NUMmodel to any FABM-enabled host model (GOTM, ROMS, SCHISM,
etc.).

---

## Directory layout

```
Fortran/fabm/
  fabm_NUMmodel.F90   — main FABM module: state variables, derivatives, sinking
  model_library.F90   — FABM model factory (registers "num/num_model")
  CMakeLists.txt      — build instructions for FABM's CMake system
  fabm.yaml.example   — annotated example configuration file
```

The NUMmodel core sources (`globals.f90`, `NUMmodel.f90`, `generalists.f90`, …)
live in `Fortran/` (one level up) and are compiled directly into the institute
library via `CMakeLists.txt`.

---

## Implementation

### Model identifier

The model is registered under the institute name **`num`** and the model name
**`num_model`**.  In `fabm.yaml` it is referenced as:

```yaml
instances:
  num_model:
    model: num/num_model
```

### State variables

| FABM flat name | Description | Units |
|---|---|---|
| `num_model_N` | Dissolved inorganic nitrogen | µgN L⁻¹ |
| `num_model_DOC` | Dissolved organic carbon | µgC L⁻¹ |
| `num_model_Si` | Dissolved silicate | µgSi L⁻¹ |
| `num_model_Gen1_1` … `Gen1_n` | Generalist biomass size classes | µgC L⁻¹ |
| `num_model_Dia1_1` … `Dia1_n` | Diatom biomass size classes | µgC L⁻¹ |
| `num_model_PCop1_1` … `PCop2_n` | Passive copepod size classes (one block per group) | µgC L⁻¹ |
| `num_model_ACop1_1` … `ACop3_n` | Active copepod size classes (one block per group) | µgC L⁻¹ |
| `num_model_POM1_1` | POM | µgC L⁻¹ |

Variable names include the per-type group counter (e.g. `PCop1_`, `PCop2_`) to
guarantee uniqueness across multiple groups of the same type.  FABM silently
aliases duplicate names, which would cause all groups of the same type to share
the same state — the naming convention here prevents that.

### Diagnostics

| FABM flat name | Description | Units |
|---|---|---|
| `num_model_ProdGross` | Gross primary production | mgC m⁻³ d⁻¹ |
| `num_model_ProdNet` | Net primary production | mgC m⁻³ d⁻¹ |
| `num_model_ProdHTL` | Production removed by higher trophic levels | mgC m⁻³ d⁻¹ |
| `num_model_BGen1` … | Total biomass per group | µgC L⁻¹ |
| `num_model_BDia1` … | | |
| `num_model_BPCop1`, `BPCop2` … | | |
| `num_model_BACop1`, `BACop2` … | | |
| `num_model_BPOM1` … | | |
| `num_model_Bpico` | Pico-plankton biomass (ESD < 2 µm) | mgC m⁻³ |
| `num_model_Bnano` | Nano-plankton biomass (2–20 µm ESD) | mgC m⁻³ |
| `num_model_Bmicro` | Micro-plankton biomass (ESD > 20 µm) | mgC m⁻³ |

### Environmental dependencies

The module reads two standard FABM variables from the host model:

- `standard_variables%temperature` (°C)
- `standard_variables%downwelling_photosynthetic_radiative_flux` (W m⁻²)

### Unit conversions

NUMmodel works in units of **days**; FABM expects **seconds**.  All rates from
`calcDerivatives` (day⁻¹) and sinking velocities (m day⁻¹) are divided by
86400 before being passed to FABM.

### Thread safety

NUMmodel uses module-level workspace arrays.  **Only one FABM instance per
executable is supported** and the `do` subroutine is **not OpenMP-safe** over
spatial grid points.

---

## Parameters (`fabm.yaml`)

| Parameter | Default | Description |
|---|---|---|
| `n_size` | 10 | Size classes per generalist / diatom group |
| `n_copepod` | 10 | Size classes per copepod group |
| `n_pom` | 1 | POM size classes |
| `n_passive` | 2 | Number of passive copepod groups |
| `n_active` | 2 | Number of active copepod groups |
| `mAdultPassive1` … | 10, 1000 µgC | Adult mass of each passive copepod group |
| `mAdultActive1` … | 10, 1000 µgC | Adult mass of each active copepod group |

See `fabm.yaml.example` for a fully annotated configuration.

---

## Building

GOTM and FABM are included as git submodules.  From the **repository root**:

```bash
git submodule update --init --recursive
```

Then configure and build (adjust NetCDF paths to match your system):

```bash
NETCDF_INC=/opt/homebrew/Cellar/netcdf-fortran/4.6.2/include
NETCDF_LIB="-L/opt/homebrew/Cellar/netcdf-fortran/4.6.2/lib -L/opt/homebrew/lib -lnetcdff -lnetcdf"

cmake -S extern/gotm -B gotm_build \
  -DFABM_BASE=$(pwd)/extern/fabm \
  -DFABM_INSTITUTES=num \
  -DFABM_NUM_BASE=$(pwd)/Fortran/fabm \
  -DNetCDF_INCLUDE_DIRS="$NETCDF_INC" \
  -DNetCDF_LIBRARIES="$NETCDF_LIB" \
  -DCMAKE_BUILD_TYPE=Release

cmake --build gotm_build -j4
```

The executable is placed at `gotm_build/gotm`.

To use NUMmodel with a different FABM-enabled host model, point that host's
CMake configuration at the same institute:

```
-DFABM_INSTITUTES=num -DFABM_NUM_BASE=/path/to/NUMmodel/Fortran/fabm
```

---

## The `input/input.yaml` path issue

NUMmodel reads all biological parameters (C:N ratio, light affinities, copepod
clearance rates, …) from a YAML file at runtime.  The path is hardcoded in
`Fortran/globals.f90`:

```fortran
character(len=19) :: inputfile='../input/input.yaml'
```

This path is **relative to the working directory** from which the host model is
launched.  It resolves correctly only when the run directory is **exactly one
level below the repository root** — as is the case for `gotm_ows_papa/`:

```
NUMmodel/              ← repository root
  input/
    input.yaml         ← biological parameters
  gotm_ows_papa/       ← run GOTM from here  →  ../input/input.yaml  ✓
  some_other_case/     ← also works if at same nesting level          ✓
  subdir/case/         ← two levels deep — path will not be found     ✗
```

**If you add a test case at a different nesting level**, either:

1. Copy or symlink `input/input.yaml` to the appropriate relative location, or
2. Change the `inputfile` variable in `globals.f90` to an absolute path or a
   configurable parameter before building.
