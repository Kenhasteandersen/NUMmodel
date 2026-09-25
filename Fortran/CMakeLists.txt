set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# The build type and optimisation flags are set in the top-level CMakeLists.txt.

# OpenMP: threads over grid cells in simulateEulerCells. Disable with -DNUM_OPENMP=OFF.
# Without OpenMP the !$omp directives are comments and everything runs serially.
option(NUM_OPENMP "Compile with OpenMP" ON)
if(NUM_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  if(OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    link_libraries(OpenMP::OpenMP_Fortran)
  else()
    message(WARNING "OpenMP not found; building without threading over grid cells.")
  endif()
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

install(TARGETS NUMmodel_matlab NUMmodel_R  DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
install(TARGETS NUMmodeltest DESTINATION ${CMAKE_CURRENT_LIST_DIR})
