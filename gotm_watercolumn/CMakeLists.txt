cmake_minimum_required(VERSION 3.12)

# ---------------------------------------------------------------------------
# Build GOTM with the NUMmodel FABM coupling.
#
# Usage:
#   mkdir build && cd build
#   cmake .. \
#       -DGOTM_BASE=<path-to-gotm-source>   \
#       -DFABM_BASE=<path-to-fabm-source>   \
#       -DCMAKE_BUILD_TYPE=Release
#   make
#
# The resulting GOTM executable will be placed in build/gotm.
#
# GOTM source: https://github.com/gotm-model/gotm
# FABM source: https://github.com/fabm-model/fabm
# ---------------------------------------------------------------------------

project(gotm_nummodel Fortran C)

# ------------------------------------------------------------------
# Resolve GOTM and FABM base directories
# ------------------------------------------------------------------
foreach(pkg GOTM FABM)
  if (NOT DEFINED ${pkg}_BASE)
    if (DEFINED ENV{${pkg}_BASE})
      set(${pkg}_BASE "$ENV{${pkg}_BASE}" CACHE PATH "Path to ${pkg} source")
    else()
      message(FATAL_ERROR
        "${pkg}_BASE is not set.\n"
        "Provide -D${pkg}_BASE=<path> on the cmake command line "
        "or set the ${pkg}_BASE environment variable.")
    endif()
  endif()
endforeach()

# ------------------------------------------------------------------
# Build FABM with NUMmodel as the biogeochemistry library
# ------------------------------------------------------------------
# Tell FABM to include our custom institute "num"
set(FABM_INSTITUTES num CACHE STRING "FABM institutes to build" FORCE)

# Point FABM to our model sources (Fortran/fabm/)
set(FABM_NUM_BASE "${CMAKE_CURRENT_SOURCE_DIR}/../Fortran/fabm"
    CACHE PATH "Path to the NUMmodel FABM sources")

add_subdirectory("${FABM_BASE}" "${CMAKE_CURRENT_BINARY_DIR}/fabm")

# ------------------------------------------------------------------
# Build GOTM with FABM enabled
# ------------------------------------------------------------------
set(GOTM_USE_FABM  ON  CACHE BOOL "Enable FABM in GOTM" FORCE)
set(GOTM_USE_NetCDF ON CACHE BOOL "Enable NetCDF output" FORCE)
set(FABM_PREFIX    "${CMAKE_CURRENT_BINARY_DIR}/fabm"
    CACHE PATH "FABM install prefix" FORCE)

add_subdirectory("${GOTM_BASE}" "${CMAKE_CURRENT_BINARY_DIR}/gotm")
