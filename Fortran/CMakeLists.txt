cmake_minimum_required(VERSION 3.10)
project(NUMmodel VERSION 0.9 LANGUAGES Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(NUMmodel SHARED
	${CMAKE_CURRENT_LIST_DIR}/input.f90
	${CMAKE_CURRENT_LIST_DIR}/globals.f90
	${CMAKE_CURRENT_LIST_DIR}/spectrum.f90
	${CMAKE_CURRENT_LIST_DIR}/generalists.f90
	${CMAKE_CURRENT_LIST_DIR}/diatoms_simple.f90
	${CMAKE_CURRENT_LIST_DIR}/diatoms.f90 
	${CMAKE_CURRENT_LIST_DIR}/generalists_csp.f90
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

