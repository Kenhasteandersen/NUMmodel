# This file is intentionally minimal.
# The correct way to build GOTM with the NUMmodel FABM coupling is to
# point cmake at the GOTM source directory and pass FABM institute flags.
#
# See README.md for the full build command, which amounts to:
#
#   mkdir build && cd build
#   cmake $GOTM_BASE                                                      \
#       -DFABM_BASE=$FABM_BASE                                            \
#       -DFABM_INSTITUTES=num                                             \
#       -DFABM_NUM_BASE=<path-to-NUMmodel>/Fortran/fabm                  \
#       -DGOTM_USE_FABM=ON                                                \
#       -DCMAKE_BUILD_TYPE=Release
#
# This file exists only to document the pattern; cmake should NOT be
# invoked with this directory as the source.
