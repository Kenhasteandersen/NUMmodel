# Major changes from version 1.0

* Changed the calculation of NPP to be based on the entire primary production (not just the assimilated one). This makes the NPP significantly higher than in v.1.0. See description in getProdNet in spectrum.f90. 

The new scheme generally increased the NPP. 

![Alt text](<NPP comparison.png>)

_Comparison of old v1.0 NPP (top), new, and difference (bottom) (units are mgC/m2/d)_