# NUMmodel
Reference implementation of the Nutrient-Unicellular-Multicellular
modelling framework.  The model is described in: Serra-Pompei et al (2020): A general size- and trait-based model of plankton communities. Progress in Oceanography (189) 102473.

The core library is written in fortran90 and is interfaced from matlab. Most of the modules are also written in matlab, so compilation of the fortran code is not needed - but it speeds up most calculations by a factor 10.

### Compiling
Use the makefile in the Fortran directory. Edit the compiler and flags to suit your operating system and compile. Compile by writing: `make lib`.

### Basic structure
There are three levels of routines: top-level, medium-level and low-level.

#### Top-level routines
* `baserunChemostat(mAdult, false)`.  Runs a chemostat version of the model and plots the output. The first argument is the adult body masses of copepods (in micro gram carbon) - send in an empty list to run only with unicellular plankton. Change the second argument to `true`to use the Fortran library.
* `baserunChemostatEuler(mAdult)`. Uses the Fortran library and simple Euler time-stepping.
* `baserunGlobal()`. Runs a global simulation with only generalists. It uses transport matrices which must be downloaded separately and placed in the directory `TMs`. Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs (choose MITgcm_2.8deg).

#### Medium-level routines
The 


