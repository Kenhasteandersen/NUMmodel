# NUMmodel
Reference implementation of the **Nutrient-Unicellular-Multicellular**
modelling framework.  The model is described in: 
* Unicellular plankton: Andersen and Visser: [From cell size and first principles to structure and function of unicellular plankton communities](https://www.biorxiv.org/content/10.1101/2022.05.16.492092v3)
* Copepods: Serra-Pompei et al (2020): [A general size- and trait-based model of plankton communities](https://www.researchgate.net/publication/346939727_A_general_size-_and_trait-based_model_of_plankton_communities "Researchgate"). 
* Application to biological pump: Serra-Pompei et al (2022) Progress in Oceanography (189) 102473: [Zooplankton trophic dynamics drive carbon export efficiency](https://www.biorxiv.org/content/10.1101/2021.03.08.434455v1 "BioRxiv").
* An [introduction to the modelling principles](https://www.youtube.com/watch?v=dHqoCqaLM8w) given by Camila Serra-Pompei. 

The core library is written in Fortran2008 and is interfaced from matlab and with a minimal frontend in  R (see http://oceanlife.dtuaqua.dk/Plankton/R).

<img width="812" alt="image" src="https://github.com/user-attachments/assets/a0fc29cc-8134-4e93-8d98-4b94dcb82f83">

https://user-images.githubusercontent.com/13268353/148120839-6bbfc0ac-69f1-445b-9b9f-b3880436bf2f.mp4

_The figure above shows a setup with only unicellular generalists run with the MIT ECCO transport matrices. The inset shows a high latitude water column at 60N, 15E, with the lower panel illustrating the trophic strategies: blue for DOC uptake (osmotrophy/bacteria), green for phototrophy, and red for phagotrophy._

### Papers using the NUM model:
* Application to the biological carbon pump: C. Serra-Pompei, B.A Ward, J. Pinti, A.W Visser, T. Kiørboe, K.H Andersen (2022): Linking plankton size spectra and community composition to carbon export and its efficiency. [Global Biogeochemical Cycles 36(5), e2021GB007275](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2021GB007275).
* Grigoratou, Maria, Camila Serra-Pompei, Adam Kemberling, and Andrew J. Pershing. Hot and hungry: [A mechanistic approach to the direct and indirect effects of marine heatwaves on plankton communities](https://assets-eu.researchsquare.com/files/rs-4194638/v1/09fb2960-27da-42ac-a1d7-3ba3ea1004fc.pdf?c=1712816635). (2024).

### Installation
The library requires a recent version of matlab (2021 or later).  Installation and compilation instructions are given in the wiki.
### Basic structure
There are three levels of routines: top-level, medium-level and low-level.  There are three model systems: an upper ocean represented as a chemostat, a water column, and a global simulation with transport matrices.
#### Top-level matlab routines
These routines run a simulation and returns the results in a `sim` structure:

* `baserunChemostat()`.  Runs a chemostat version of the model and plots the output.
* `baserunWatercolumn`. Runs a watercolumn extracted from a transport matrix at a specific location.
* `baserunGlobal()`. Runs a global simulation with only unicellular generalists. 

All units are in micro gC/l (or micro gN/l for nutrient concentration). Units of light are micro mol photons per m2 per second.

#### Medium-level matlab routines
The routines operates with two basic structures: a *parameter* structure and a *simulation* structure. The parameter structure contains all parameters needed for a simulation. The simulation structure contains all the output, which can be used for analysis or for plotting.  See more details on the wiki page.

*Parameters* are set with two calls: one to setup the size spectra to simulate and one to add the parameters for the simulation (chemostat or global). The size spectra are setup with a call to `setupXX` where XX represent the setup, e.g., `setupNUMmodel` or `setupGeneralistsOnly` or (the latter onlu includes unicellular generalists).  Parameters for the simulation are subsequently set with a call to `parametersChemostat`, `parametersWatercolumn`, or `parametersGlobal`. For example: `p = parametersChemostat( setupGeneralistsOnly() );` (see the wiki for a description of the parameter structure).

*Simulations* are performed with calls to a simulation routine: `sim = simulationChemostat(p)`, `sim = simulationWatercolumn(p, latitude, longitude)`, or `sim = simulationGlobal(p)`, where `p` is the parameter structure (see the wiki for a description of the simulation structure).

*Plots* are made with calls to the plot routines. `plotSimulation(sim)` makes a series of basic plots of a simulation.
