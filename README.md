# NUMmodel
Reference implementation of the **Nutrient-Unicellular-Multicellular**
modelling framework.  The model is described in: Serra-Pompei et al (2020): [A general size- and trait-based model of plankton communities](https://www.researchgate.net/publication/346939727_A_general_size-_and_trait-based_model_of_plankton_communities "Available on Researchgate"). Progress in Oceanography (189) 102473.

:fire: **This is alpha-software in development. Things may be broken and results may not be correct** :fire:

The core library is written in fortran90 and is interfaced from matlab. Most of the modules are also written in matlab so compilation of the fortran code is not needed - but it speeds up most calculations by a factor 10.

### Compiling
Use the makefile in the Fortran directory. Edit the compiler and flags to suit your operating system and compile. Compile by writing: `make lib`.

### Basic structure
There are three levels of routines: top-level, medium-level and low-level.  There are two model systems: an upper ocean represented as a chemostat and a global simulation with transport matrices.  See `exampleGeneralists` for some basic runs of the chemostat model.

#### Top-level routines
* `baserunChemostat(mAdult, false)`.  Runs a chemostat version of the model and plots the output. The first argument is the adult body masses of copepods (in micro gram carbon) - send in an empty list to run only with unicellular plankton. Change the second argument to `true`to use the Fortran library.
* `baserunChemostatEuler(mAdult)`. Uses the Fortran library and simple Euler time-stepping.
* `baserunGlobal()`. Runs a global simulation with only generalists. It uses transport matrices which must be downloaded separately and placed in the directory `TMs`. Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs (choose MITgcm_2.8deg).

#### Medium-level routines
The routines operates with two basic structures: a *parameter* structure and a *simulation* structure. The parameter structure contains all parameters needed for a simulation. The simulation structure contains all the output, which can be used for analysis or for plotting.

*Parameters* are set first with a call to `parameters(mAdult)` where `mAdult` again is vector of copepod adult masses.  Additional parameters are subsequently set with a call to `parametersChemostat`or `parametersGlobal`. For example: `p = parametersChemostat( parameters ([]) );`.

*Simulations* are performed with calls to a simulation routine: `sim = simulationChemostat(p)` or `sim = simulationGlobal(p)`, where `p` is a parameter structure.

*Plots* are made with calls to the plot routines, e.g. `plotChemostat(sim)` for a chemostat run.


