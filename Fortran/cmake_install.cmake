# Install script for directory: /home/trine/Documents/GitHub/NUMmodel/Fortran

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib" TYPE SHARED_LIBRARY FILES "/home/trine/Documents/GitHub/NUMmodel/Fortran/libNUMmodel_matlab.so")
  if(EXISTS "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_matlab.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib" TYPE SHARED_LIBRARY FILES "/home/trine/Documents/GitHub/NUMmodel/Fortran/libNUMmodel_R.so")
  if(EXISTS "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/../lib/libNUMmodel_R.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/trine/Documents/GitHub/NUMmodel/Fortran" TYPE EXECUTABLE FILES "/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest")
  if(EXISTS "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/trine/Documents/GitHub/NUMmodel/Fortran/NUMmodeltest")
    endif()
  endif()
endif()

