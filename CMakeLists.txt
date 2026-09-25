cmake_minimum_required(VERSION 3.10)

# Avoid CMake picking up an unlicensed NAG Fortran compiler (nagfor) ahead of
# gfortran when both are installed; CMake's default search prefers "paid"
# compilers like nagfor over gfortran.
if(NOT CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

project(NUMmodel VERSION 0.9 LANGUAGES Fortran)

# Build optimised by default. Without an explicit -O flag gfortran compiles at
# -O0, which makes the library about 4.5 times slower. -O3 measured no faster
# than -O2 for this code, so Release is pinned to -O2.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -DNDEBUG")

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Fortran)