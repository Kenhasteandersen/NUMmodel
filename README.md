# NUMmodel
Reference implementation of the **Nutrient-Unicellular-Multicellular**
modelling framework.  The model is described in: Serra-Pompei et al (2020): [A general size- and trait-based model of plankton communities](https://www.researchgate.net/publication/346939727_A_general_size-_and_trait-based_model_of_plankton_communities "Researchgate"). Progress in Oceanography (189) 102473 and Serra-Pompei et al: [Zooplankton trophic dynamics drive carbon export efficiency](https://www.biorxiv.org/content/10.1101/2021.03.08.434455v1 "BioRxiv").

:fire: **This is alpha-software in development. Things may be broken and results may not be correct** :fire:

The core library is written in fortran90 and is interfaced from matlab or R. Requires a fairly recent version of matlab to run.

### Compiling
Use the makefile in the Fortran directory. Edit the compiler and flags in the makefile to suit your operating system and compile by writing: `make`.

### Basic structure
There are three levels of routines: top-level, medium-level and low-level.  There are two model systems: an upper ocean represented as a chemostat and a global simulation with transport matrices.
#### Top-level matlab routines
These routines run a simulation and returns the results in a `sim` structure. All units are in micro gC/l (or micro gN/l). Units of light are micro mol photons per m2 per second:

* `baserunChemostat(mAdult)`.  Runs a chemostat version of the model and plots the output. The argument is the adult body masses of copepods (in micro gram carbon) - send in an empty list to run only with unicellular plankton.
* `baserunChemostatEuler(mAdult)`. Uses the Fortran library and simple Euler time-stepping.
* `baserunGlobal()`. Runs a global simulation with only generalists. It uses transport matrices which must be downloaded separately and placed in the directory `TMs`. Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs (choose MITgcm_2.8deg).
* `baserunWatercolumn`. Run a watercolumn extracted from a transport matrix.

#### Medium-level matlab routines
The routines operates with two basic structures: a *parameter* structure and a *simulation* structure. The parameter structure contains all parameters needed for a simulation. The simulation structure contains all the output, which can be used for analysis or for plotting.

*Parameters* are set with two calls: one to setup the size spectra to simulate and one to add the parameters for the simulation (chemostat or global). The size spectra are setup with a call to `setupXX` where XX represent the setup, e.g., `setupGeneralistsOnly` or `setupGeneric(mAdult)` (the latter includes copepods where `mAdult` is vector of copepod adult masses).  Parameters for the simulation are subsequently set with a call to `parametersChemostat`or `parametersGlobal`. For example: `p = parametersChemostat( setupGeneralistsOnly() );`.

*Simulations* are performed with calls to a simulation routine: `sim = simulationChemostat(p)`, `sim = simulationWatercolumn(p, latitude, longitude)`, or `sim = simulationGlobal(p)`, where `p` is a parameter structure.

*Plots* are made with calls to the plot routines. A series of basic plots are made by a call to `plotSimulation(sim)`.
