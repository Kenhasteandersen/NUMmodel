cmake_minimum_required(VERSION 3.12)

# ---------------------------------------------------------------------------
# FABM coupling for NUMmodel
#
# Build instructions:
#   1. Clone FABM: git clone https://github.com/fabm-model/fabm.git
#   2. From this directory:
#        mkdir build && cd build
#        cmake .. \
#          -DFABM_BASE=<path-to-fabm> \
#          -DFABM_HOST=<host-driver>   # e.g. gotm, roms, 0d
#        make
#
# The CMake variable FABM_BASE must point to the root of the FABM source tree.
# Alternatively, set the environment variable FABM_BASE before running CMake.
# ---------------------------------------------------------------------------

project(fabm_num Fortran)

# Locate FABM base directory (from CMake variable or environment)
if (NOT DEFINED FABM_BASE)
  if (DEFINED ENV{FABM_BASE})
    set(FABM_BASE "$ENV{FABM_BASE}" CACHE PATH "Path to FABM source")
  else()
    message(FATAL_ERROR
      "FABM_BASE is not set. "
      "Please provide -DFABM_BASE=<path-to-fabm> on the cmake command line, "
      "or set the FABM_BASE environment variable.")
  endif()
endif()

# Include FABM's CMake infrastructure
add_subdirectory("${FABM_BASE}" "${CMAKE_CURRENT_BINARY_DIR}/fabm_base" EXCLUDE_FROM_ALL)

# ---------------------------------------------------------------------------
# NUMmodel core source files (relative to this CMakeLists.txt)
# ---------------------------------------------------------------------------
set(NUMMODEL_SOURCES
    ../globals.f90
    ../read_input.f90
    ../spectrum.f90
    ../generalists_simple.f90
    ../generalists.f90
    ../diatoms_simple.f90
    ../diatoms.f90
    ../copepods.f90
    ../POM.f90
    ../NUMmodel.f90
)

# ---------------------------------------------------------------------------
# FABM coupling sources
# ---------------------------------------------------------------------------
set(FABM_NUM_SOURCES
    fabm_NUMmodel.F90
    model_library.F90
)

# Use FABM's helper to create the model library
fabm_add_library(fabm_num
    ${NUMMODEL_SOURCES}
    ${FABM_NUM_SOURCES}
)

# The NUMmodel input.yaml is read at runtime from ../input/input.yaml
# (relative to the working directory). No extra install step needed here,
# but the user must ensure the input/ directory is accessible at runtime.
