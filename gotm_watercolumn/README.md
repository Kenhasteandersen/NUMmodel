# NUMmodel GOTM water column example

A 100 m North Sea water column (59°N, 0°E) running for one year (2000),
using GOTM for the physical oceanography and FABM to call NUMmodel for
the biogeochemistry.

## Files

| File | Description |
|---|---|
| `gotm.yaml` | GOTM main configuration |
| `fabm.yaml` | FABM / NUMmodel configuration |
| `init_temp.dat` | Initial winter temperature profile |
| `generate_forcing.py` | Script that generates `meteo.dat` |
| `CMakeLists.txt` | Build system for GOTM+FABM+NUMmodel |

## Prerequisites

- CMake ≥ 3.12
- gfortran ≥ 9 (or ifort ≥ 19)
- NetCDF with Fortran bindings (`libnetcdff`)
- Python 3 with NumPy (for `generate_forcing.py`)

Clone the required external models:

```bash
# GOTM
git clone https://github.com/gotm-model/gotm.git  ~/src/gotm

# FABM
git clone https://github.com/fabm-model/fabm.git  ~/src/fabm
```

## Build

### 1. Register the NUMmodel FABM institute

FABM discovers model libraries via a `${institute}_model_library.F90` convention.
The FABM coupling sources are in `Fortran/fabm/` (one directory above this one).

FABM needs the NUMmodel institute sources on its include path.
The `CMakeLists.txt` in this directory handles this automatically by setting
`FABM_INSTITUTES=num` and pointing `FABM_NUM_BASE` to `../Fortran/fabm/`.

### 2. Configure and compile

```bash
cd NUMmodel/gotm_watercolumn
mkdir build && cd build

cmake ..                              \
  -DGOTM_BASE=~/src/gotm             \
  -DFABM_BASE=~/src/fabm             \
  -DCMAKE_BUILD_TYPE=Release

make -j$(nproc)
```

The GOTM executable is at `build/gotm/gotm` (or `build/gotm/src/gotm`
depending on the GOTM version).

### 3. Verify

```bash
./build/gotm/gotm --version
```

## Run

### 1. Generate meteorological forcing

```bash
cd NUMmodel/gotm_watercolumn
python generate_forcing.py
# Creates meteo.dat (365 daily records)
```

### 2. Run GOTM

GOTM must be run from the `gotm_watercolumn/` directory so that
NUMmodel can find its parameter file at the relative path `../input/input.yaml`.

```bash
cd NUMmodel/gotm_watercolumn
./build/gotm/gotm
```

GOTM reads `gotm.yaml` and `fabm.yaml` from the current directory.
Output is written to `NUMmodel_watercolumn.nc` (NetCDF).

## Output

The NetCDF file contains daily snapshots of all GOTM physical variables
plus all NUMmodel biogeochemical state variables and diagnostics.

Key variables to examine first:

| GOTM variable | Description | Units |
|---|---|---|
| `temp` | Water temperature | °C |
| `num_model_N` | Dissolved inorganic nitrogen | µgN/L |
| `num_model_Si` | Dissolved silicate | µgSi/L |
| `num_model_Gen1..10` | Generalist biomass size classes | µgC/L |
| `num_model_Dia1..10` | Diatom biomass size classes | µgC/L |
| `num_model_PCop1..10` | Passive copepod size classes | µgC/L |
| `num_model_ACop1..10` | Active copepod size classes | µgC/L |
| `num_model_POM1` | Particulate organic matter | µgC/L |
| `num_model_ProdGross` | Gross primary production | mgC/d/m³ |
| `num_model_Bpico` | Pico-plankton biomass | mgC/m³ |

### Quick Python plot

```python
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

ds = xr.open_dataset("NUMmodel_watercolumn.nc")

# Total generalist biomass (sum over size classes)
gen_vars = [v for v in ds if v.startswith("num_model_Gen")]
B_gen = sum(ds[v] for v in gen_vars)

# Total diatom biomass
dia_vars = [v for v in ds if v.startswith("num_model_Dia")]
B_dia = sum(ds[v] for v in dia_vars)

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# Temperature
ds["temp"].plot(ax=axes[0, 0], cmap="RdYlBu_r")
axes[0, 0].set_title("Temperature (°C)")

# Dissolved nitrogen
ds["num_model_N"].plot(ax=axes[0, 1], cmap="YlOrRd")
axes[0, 1].set_title("DIN (µgN/L)")

# Generalist biomass
B_gen.plot(ax=axes[1, 0], cmap="Greens")
axes[1, 0].set_title("Generalist biomass (µgC/L)")

# Diatom biomass
B_dia.plot(ax=axes[1, 1], cmap="Blues")
axes[1, 1].set_title("Diatom biomass (µgC/L)")

plt.tight_layout()
plt.savefig("NUMmodel_watercolumn.png", dpi=150)
plt.show()
```

## Tuning

All ecological parameters (affinity constants, predator–prey ratios, etc.)
are read from `../input/input.yaml`.  Biological initial conditions and
the copepod adult-mass configuration are set in `fabm.yaml`.

The number of size classes and copepod groups can be changed in `fabm.yaml`:
```yaml
    parameters:
      n_size:    10   # size classes per generalist/diatom group
      n_copepod: 10   # size classes per copepod group
      n_pom:      1   # POM size classes
      n_passive:  2   # passive copepod groups
      n_active:   2   # active copepod groups
```

Reducing `n_size` and `n_copepod` to 5 significantly speeds up the run.

## Known limitations

- NUMmodel uses module-level global arrays, so only **one FABM instance**
  per executable is supported; do not declare multiple `num/num_model` instances.
- Thread-safe OpenMP parallelism over spatial points is not guaranteed.
  Run GOTM with `OMP_NUM_THREADS=1` if you encounter race conditions.
- The NUMmodel parameter file path is hardcoded as `../input/input.yaml`
  (relative to the GOTM working directory).  Always run GOTM from the
  `gotm_watercolumn/` directory, or adjust the path in `Fortran/globals.f90`.
