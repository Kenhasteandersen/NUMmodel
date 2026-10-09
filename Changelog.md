# Major changes from version 1.0

* Changed the calculation of NPP to be based on the entire primary production (not just the assimilated one).  This makes the NPP significantly higher than in v.1.0. See description in getProdNet in spectrum.f90.  The new scheme generally increased the NPP; see below.

![Alt text](<NPP comparison.png>)

_Comparison of old v1.0 NPP (top), new, and difference (bottom) (units are mgC/m2/d)_

* Re-calibrated to fit with the new NPP scheme: kw=0.10/m, u=20 m/d, mortHTL=0.006.

* Several fixed resulted in factor 3 or more speed increases.

* Added support for UVic transport matrix

* Fixed loss of N in the transport matrices. The global model now uses the raw transport matrices and removes the negative concentrations they create after each transport step, by taking the missing mass from the nearest boxes above and below in the same water column (function_clip_negative; p.bClipNegative = true). The previous conversion of the matrices to non-negative ones conserved mass but did not keep a uniform field uniform, and it moved tracer from small to large boxes (downwards). Mass that a water column cannot cover is added and reported in sim.Nclip and sim.clipStats. Set p.bClipNegative = false for the old behaviour.

# Minor fixes and changes:

* Fixed simulateChemostatEuler to give similar results as simulateChemostat

* Fixed a bug in the predictor-corrector scheme for DOC

* Fixed the calculation of N balance in the water column, chemostat, and global simulations. Nothing in the model has been changed; it was only the balance calculations which were incorrect.

* Fixed the conservation of mass in the explict part of the transport matrix.

* Moved the input file to a .yaml file
