# Minimal CMakeLists.txt for ALPS embedded into QCMaquis
# (c) 2020 Leon Freitag

## Original ALPS CMakeLists.txt header:
#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


######################################################################
# CMAKE_MODULE_PATH
######################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

######################################################################
# Version information
######################################################################
set(ALPS_YEAR 2016)
set(ALPS_VERSION_MAJOR 2)
set(ALPS_VERSION_MINOR 3)
set(ALPS_VERSION_PATCH 0)
set(ALPS_VERSION_BUILD "")

IF(NOT ALPS_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn)
  FIND_PACKAGE(Subversion)
  IF(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} ALPS)
    set(ALPS_VERSION_BUILD "r${ALPS_WC_REVISION}")
  ENDIF(Subversion_FOUND)
ENDIF(NOT ALPS_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn)

if(ALPS_VERSION_BUILD)
  set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}-${ALPS_VERSION_BUILD}")
else(ALPS_VERSION_BUILD)
  set(ALPS_VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}")
endif(ALPS_VERSION_BUILD)
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")
MESSAGE(STATUS "ALPS version: ${ALPS_VERSION}")

set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)

SET(ALPS_HAVE_MPI 0)
SET(ALPS_HAVE_BOOST_MPI 0)
set(BUILD_BOOST_MPI FALSE)
# endif(ALPS_ENABLE_MPI)

# Boost include, should have been found by find_boost in parent CMakeLists.txt
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
link_libraries(${Boost_LIBRARIES})

IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  add_definitions(${HDF5_DEFINTIONS})
  SET(ALPS_HAVE_HDF5 1)
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "HDF5 >= 1.8.2 or higher has not been found")
ENDIF(HDF5_FOUND)

IF(HDF5_IS_PARALLEL)
  SET(ALPS_HAVE_HDF5_PARALLEL 1)
  MESSAGE(WARNING "parallel(MPI) hdf5 is detected. We will compile but ALPS does not use parallel HDF5. The standard version is preferred.")
  IF(NOT MPI_FOUND)
    MESSAGE(FATAL_ERROR "parallel(MPI) hdf5 needs MPI. Enable MPI or install serial HDF5 libraries.")
  ENDIF(NOT MPI_FOUND)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF(HDF5_IS_PARALLEL)


######################################################################
# configuration files
######################################################################

# This should be covered by configure_file, but for some reason it isn't
# probably because the configure file is included at a wrong place
# so let's set it again explicitly
#
if(LAPACK_64_BIT)
   add_definitions(-DBIND_FORTRAN_INTEGER_8)
endif()
if (BIND_FORTRAN_LOWERCASE)
  add_definitions(-DBIND_FORTRAN_LOWERCASE)
endif (BIND_FORTRAN_LOWERCASE)

set(ALPS_XML_PATH "lib")
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/${ALPS_XML_PATH}/xml")

configure_file(src/alps/config.h.in ${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src/alps/config.h)
configure_file(src/alps/version.h.in ${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src/alps/version.h)

# order the include directories so that directories which are in the
# source or build tree always come before directories outside the project.
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

######################################################################
# include directories
######################################################################
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src)
include_directories(src)

######################################################################
# RPATH setting
######################################################################
if(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_MACOSX_RPATH 1)
endif(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))

## Disabled tests
# if (ALPS_BUILD_TESTS AND NOT ALPS_BUILD_SOURCE)
#   include(ALPSBuildname)
#   enable_testing()
#   include(CTest)
#   add_subdirectory(test)
# endif (ALPS_BUILD_TESTS AND NOT ALPS_BUILD_SOURCE)
# add_subdirectory(config)
add_subdirectory(src/alps)

