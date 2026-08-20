# Major changes from version 1.0

* Changed the calculation of NPP to include the primary production that is lost to
DOC, but then quickly taken up again. This makes the NPP significantly higher than in v.1.0. See description in getProdNet in spectrum.f90. 

The new scheme generally increased the NPP. 

![Alt text](NPP\ comparison.png)
_comparison of old v1.0 NPP (top), new, and difference (bottom)