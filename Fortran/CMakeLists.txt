set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Specify default build type for single-build-type systems (not VS)
#if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
#   set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
#   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
#endif()
set(CMAKE_FORTRAN_FLAGS="-O2")
set(RELEASE_FLAGS "-O2 -DNDEBUG -DBOOST_DISABLE_ASSERTS")
set(DEBUG_FLAGS "-g3")

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
