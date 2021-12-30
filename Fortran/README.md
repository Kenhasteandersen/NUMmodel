Files for the Fortran library. Requires a F2008 compatible compiler. Tested with gfortran and intel fortran.

Files:
 - `NUMmodel.f90`. The core library
 - `NUMmodel_wrap_colmajor.f90`. Wrapper for matlab
 - `NUMmodel_wrap_R.f90`. Wrapper for R. Does not implement all function in the library
 - `globals.f90`. Global definitions
 - `spectrum.f90`. The abstract class for all spectrum groups. Contains definitions of the parent `typeSpectrum` class and the two derived classes `spectrumUnicellular` and `spectrumMulticellular`.

 The library current implements the following plankton groups:
- `Generalists.f90`. Mixotrophs (bacteria, flaggelates, dinoflaggelates and ciliates). Tested and working.
- `Copepods.f90`. Copepods defined by their adult size. Working but not thoroughly tested.
- `Diatoms.f90`. Diatoms. Not working.
- `Diatoms_simple.f90`. Very simple version of diatoms. Working but not thoroughly tested.
- `Generalists_csp.f90`. Implementation of the simpler version of generalists from Serra-Pompei (2020). Not working.