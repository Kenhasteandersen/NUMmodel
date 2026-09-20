cmake_minimum_required(VERSION 3.10)

# Avoid CMake picking up an unlicensed NAG Fortran compiler (nagfor) ahead of
# gfortran when both are installed; CMake's default search prefers "paid"
# compilers like nagfor over gfortran.
if(NOT CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

project(NUMmodel VERSION 0.9 LANGUAGES Fortran)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Fortran)