set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Specify default build type for single-build-type systems (not VS)
#if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#   set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
#   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
#endif()
set(CMAKE_FORTRAN_FLAGS="-O2")
set(RELEASE_FLAGS "-O2 -DNDEBUG -DBOOST_DISABLE_ASSERTS")
set(DEBUG_FLAGS "-g3 -gline")

# OpenMP: threads over grid cells in simulateEulerCells. Disable with -DNUM_OPENMP=OFF
option(NUM_OPENMP "Compile with OpenMP" ON)
if(NUM_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  if(OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    link_libraries(OpenMP::OpenMP_Fortran)
  endif()
endif()

# GPU offload of the prototype kernel in NUMmodel_offload.f90. Examples:
#   gfortran:  -DNUM_OFFLOAD=ON -DNUM_OFFLOAD_FLAGS="-foffload=nvptx-none"   (or amdgcn-amdhsa)
#   nvfortran: -DNUM_OFFLOAD=ON -DNUM_OFFLOAD_FLAGS="-mp=gpu -gpu=cc80"
#   amdflang:  -DNUM_OFFLOAD=ON -DNUM_OFFLOAD_FLAGS="--offload-arch=gfx90a"
#   ifx:       -DNUM_OFFLOAD=ON -DNUM_OFFLOAD_FLAGS="-fopenmp-targets=spir64"
# With NUM_OFFLOAD=OFF target regions run on the host.
option(NUM_OFFLOAD "Compile OpenMP target regions for a GPU" OFF)
set(NUM_OFFLOAD_FLAGS "" CACHE STRING "Compiler flags for OpenMP offload")
if(NUM_OFFLOAD)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NUM_OFFLOAD_FLAGS}")
  link_libraries(${NUM_OFFLOAD_FLAGS})
elseif(NUM_OPENMP AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # gcc builds with offload support otherwise try to add device code
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -foffload=disable")
endif()

add_library(NUMmodel SHARED
	${CMAKE_CURRENT_LIST_DIR}/read_input.f90
	${CMAKE_CURRENT_LIST_DIR}/globals.f90
	${CMAKE_CURRENT_LIST_DIR}/spectrum.f90
    ${CMAKE_CURRENT_LIST_DIR}/generalists_simple.f90
	${CMAKE_CURRENT_LIST_DIR}/generalists.f90
	${CMAKE_CURRENT_LIST_DIR}/diatoms_simple.f90
	${CMAKE_CURRENT_LIST_DIR}/diatoms.f90 
	${CMAKE_CURRENT_LIST_DIR}/copepods.f90
	${CMAKE_CURRENT_LIST_DIR}/POM.f90
	${CMAKE_CURRENT_LIST_DIR}/NUMmodel.f90
	${CMAKE_CURRENT_LIST_DIR}/NUMmodel_offload.f90
)

add_library(NUMmodel_matlab SHARED
	$<TARGET_OBJECTS:NUMmodel>
	${CMAKE_CURRENT_LIST_DIR}/NUMmodel_wrap_colmajor.f90
)

add_library(NUMmodel_R SHARED
    $<TARGET_OBJECTS:NUMmodel>
    ${CMAKE_CURRENT_LIST_DIR}/NUMmodel_wrap_R.f90
)

add_executable(NUMmodeltest 
	$<TARGET_OBJECTS:NUMmodel>
	${CMAKE_CURRENT_LIST_DIR}/NUMmodeltest.f90
)

add_executable(NUMmodeltest_cells
	$<TARGET_OBJECTS:NUMmodel>
	${CMAKE_CURRENT_LIST_DIR}/NUMmodeltest_cells.f90
)

install(TARGETS NUMmodel_matlab NUMmodel_R  DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
install(TARGETS NUMmodeltest DESTINATION ${CMAKE_CURRENT_LIST_DIR})
