# NUMmodel / GOTM — OWS Papa test case

Water-column simulation of the North Pacific at Ocean Weather Station Papa
(50°N, 145°W) using the [GOTM](https://gotm.net) turbulence model coupled to
the NUMmodel size-structured plankton model via
[FABM](https://github.com/fabm-model/fabm).

Observed meteorological and hydrographic forcing covers **2011–2020** (9 years).

---

## Prerequisites

| Tool | Tested version |
|------|---------------|
| gfortran | 15 (Homebrew) |
| CMake | ≥ 3.20 |
| NetCDF-Fortran | 4.6.2 |
| Python 3 | 3.14 (Homebrew) |
| Python packages | `netCDF4`, `matplotlib`, `numpy` |

On macOS with Homebrew:
```bash
brew install gcc cmake netcdf netcdf-fortran
pip3 install netCDF4 matplotlib numpy
```

---

## 1. Clone the repository with submodules

GOTM and FABM are included as git submodules and are **not** bundled in this
repository — only the commit hashes are stored here.

```bash
git clone https://github.com/Kenhasteandersen/NUMmodel.git
cd NUMmodel
git submodule update --init --recursive
```

`--recursive` is needed because GOTM itself has nested submodules
(GSW-Fortran, flexout, CVMix, etc.).

After this step the directory layout is:
```
NUMmodel/
  extern/
    gotm/       ← GOTM source (https://github.com/gotm-model/code.git)
    fabm/       ← FABM source (https://github.com/fabm-model/fabm.git)
  Fortran/
    fabm/       ← NUMmodel FABM coupling (this repo)
  gotm_ows_papa/   ← this directory
  gotm_build/      ← created in step 2 (not tracked by git)
```

---

## 2. Build GOTM with the NUMmodel FABM institute

Run the following from the **repository root** (`NUMmodel/`):

```bash
cd NUMmodel   # repository root

# Adjust NetCDF paths to match your system
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

On success the executable is at `gotm_build/gotm`.

### Finding your NetCDF paths

```bash
# Include directory
nc-config --includedir       # or: nf-config --includedir

# Library directory
nc-config --libdir
nf-config --libdir
```

---

## 3. Run the test case

```bash
cd gotm_ows_papa
python3 run_and_plot.py
```

The script:
1. Runs `gotm_build/gotm` in the `gotm_ows_papa/` directory
2. Reads the NetCDF output `ows_papa_NUMmodel.nc`
3. Saves three plots:

| File | Contents |
|------|----------|
| `all_variables.png` | Depth–time colour panels for all 19 state variables |
| `integrated_variables.png` | Depth-integrated time series (nutrients + biomass groups) |
| `diffusivity.png` | Turbulent diffusivity $K_h$ on a log colour scale |

If your GOTM executable is in a different location:
```bash
python3 run_and_plot.py --gotm /path/to/gotm
```

---

## Configuration files

| File | Purpose |
|------|---------|
| `gotm.yaml` | GOTM physics configuration (grid, turbulence, surface fluxes, output) |
| `fabm.yaml` | NUMmodel biology configuration (group sizes, copepod adult masses, initial conditions) |

### Physical setup (`gotm.yaml`)

- **Location:** 50°N, 145°W, 150 m depth, 150 vertical layers
- **Period:** 2011-03-21 – 2020-03-21 (9 years)
- **Turbulence:** k-ε second-order closure with Canuto-A stability functions and
  Large et al. (1994) internal wave mixing
- **Forcing:** observed hourly time series of temperature profiles, salinity
  profiles (with daily relaxation), wind stress, shortwave radiation, longwave
  radiation, and turbulent heat fluxes from the standard OWS Papa dataset

### Biological setup (`fabm.yaml`)

The NUMmodel groups match the MATLAB `setupNUMmodel` defaults:

| Group | Setting |
|-------|---------|
| Generalists | 1 group × 10 size classes |
| Diatoms | 1 group × 10 size classes |
| Passive copepods | 2 groups (adult mass 0.2 and 5 µgC), 6 size classes each |
| Active copepods | 3 groups (adult mass 1, 31.6 and 1000 µgC), 6 size classes each |
| POM | 1 class |

Initial nutrient conditions: N = 150 µgN/L, Si = 200 µgSi/L, DOC = 0.

---

## Forcing data

The forcing files (`tprof_papa_hourly.dat`, `sprof_papa_hourly.dat`,
`heat_flux_papa.dat`, `momentum_flux_papa.dat`, `swr_papa.dat`, `lwr.dat`,
`airt.dat`, `hum.dat`, `airp.dat`, `u10.dat`, `sst_hourly.dat`,
`sss_hourly.dat`) are derived from the standard GOTM OWS Papa test case
available at <https://github.com/gotm-model/cases>.
