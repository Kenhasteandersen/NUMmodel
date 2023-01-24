# NUMmodel
Reference implementation of the **Nutrient-Unicellular-Multicellular**
modelling framework.  The model is described in: 
* Unicellular plankton: Andersen and Visser: [From cell size and first principles to structure and function of unicellular plankton communities](https://www.biorxiv.org/content/10.1101/2022.05.16.492092v3)
* Copepods: Serra-Pompei et al (2020): [A general size- and trait-based model of plankton communities](https://www.researchgate.net/publication/346939727_A_general_size-_and_trait-based_model_of_plankton_communities "Researchgate"). 
* Application to biological pump: Serra-Pompei et al (2022) Progress in Oceanography (189) 102473: [Zooplankton trophic dynamics drive carbon export efficiency](https://www.biorxiv.org/content/10.1101/2021.03.08.434455v1 "BioRxiv").

Try the [online simulator of the unicellular model](http://oceanlife.dtuaqua.dk/Plankton/R/).

:fire: **This is alpha-software in development. Things may be broken and results may not be correct** :fire:



https://user-images.githubusercontent.com/13268353/148120839-6bbfc0ac-69f1-445b-9b9f-b3880436bf2f.mp4




The core library is written in Fortran2008 and is interfaced from matlab or R.

### Installation
The library requires a recent version of matlab (2021 or later). On windows it requires the Matlab MEX module to be installed (Home -> Add-ons -> Get Add-ons -> MATLAB Support for MinGW-w64 C/C++ Compiler); on mac it requires Xcode to be installed. To run global simulation it further requires that the mapping toolbox is installed.  Compiled versions of the library is available for windows (64 bit), linux and osx.  Compiling the library requires a Fortran compiler, e.g., gfortran.  Use the makefile in the Fortran directory. Edit the compiler and flags in the makefile to suit your operating system and compile by writing: `make`.

### Basic structure
There are three levels of routines: top-level, medium-level and low-level.  There are two model systems: an upper ocean represented as a chemostat and a global simulation with transport matrices.
#### Top-level matlab routines
These routines run a simulation and returns the results in a `sim` structure:

* `baserunChemostat(mAdult)`.  Runs a chemostat version of the model and plots the output. The argument is the adult body masses of copepods (in micro gram carbon) - send in an empty list to run only with unicellular plankton.
* `baserunChemostatEuler(mAdult)`. Uses the Fortran library and simple Euler time-stepping.
* `baserunGlobal()`. Runs a global simulation with only generalists. It uses transport matrices which must be downloaded separately and placed in the directory `TMs`. Transport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs (choose MITgcm_2.8deg).
* `baserunWatercolumn`. Runs a watercolumn extracted from a transport matrix at a specific location.

All units are in micro gC/l (or micro gN/l for nutrient concentration). Units of light are micro mol photons per m2 per second.

#### Medium-level matlab routines
The routines operates with two basic structures: a *parameter* structure and a *simulation* structure. The parameter structure contains all parameters needed for a simulation. The simulation structure contains all the output, which can be used for analysis or for plotting.

*Parameters* are set with two calls: one to setup the size spectra to simulate and one to add the parameters for the simulation (chemostat or global). The size spectra are setup with a call to `setupXX` where XX represent the setup, e.g., `setupGeneralistsOnly` or `setupGeneric(mAdult)` (the latter includes copepods where `mAdult` is vector of copepod adult masses).  Parameters for the simulation are subsequently set with a call to `parametersChemostat`or `parametersGlobal`. For example: `p = parametersChemostat( setupGeneralistsOnly() );`.

*Simulations* are performed with calls to a simulation routine: `sim = simulationChemostat(p)`, `sim = simulationWatercolumn(p, latitude, longitude)`, or `sim = simulationGlobal(p)`, where `p` is a parameter structure.

*Plots* are made with calls to the plot routines. A series of basic plots are made by a call to `plotSimulation(sim)`.
