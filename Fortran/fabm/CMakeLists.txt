# FABM institute library for NUMmodel ("num")
#
# This CMakeLists.txt is invoked by FABM's build system when
# FABM_INSTITUTES=num and FABM_NUM_BASE points to this directory.
# It must create an OBJECT library named fabm_models_num.

add_library(fabm_models_num OBJECT
    # FABM coupling
    fabm_NUMmodel.F90
    model_library.F90
    # NUMmodel core (relative paths to Fortran/ directory)
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

target_link_libraries(fabm_models_num PRIVATE fabm_base)
