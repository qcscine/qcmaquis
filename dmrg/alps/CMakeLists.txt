# Minimal CMakeLists.txt for ALPS embedded into QCMaquis (c) 2020 Leon Freitag

# Original ALPS CMakeLists.txt header: Copyright Synge Todo and Matthias Troyer
# 2009 - 2010. Distributed under the Boost Software License, Version 1.0. (See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# ##############################################################################
# CMAKE_MODULE_PATH
# ##############################################################################
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# ##############################################################################
# Version information
# ##############################################################################
set(ALPS_YEAR 2016)
set(ALPS_VERSION_MAJOR 2)
set(ALPS_VERSION_MINOR 3)
set(ALPS_VERSION_PATCH 0)
set(ALPS_VERSION_BUILD "")
set(ALPS_VERSION
    "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}")
set(ALPS_VERSION_STRING "ALPS Libraries version ${ALPS_VERSION}")
set(ALPS_CONFIG_HOST unknown)
set(ALPS_CONFIG_USER unknown)

set(ALPS_HAVE_MPI 0)
set(ALPS_HAVE_BOOST_MPI 0)
set(BUILD_BOOST_MPI FALSE)

# Boost include, should have been found by find_boost in parent CMakeLists.txt
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
link_libraries(${Boost_LIBRARIES})

if(HDF5_FOUND)
  set(HAVE_LIBHDF5 1)
  include_directories(${HDF5_INCLUDE_DIR})
  add_definitions(${HDF5_DEFINTIONS})
  set(ALPS_HAVE_HDF5 1)
else(HDF5_FOUND)
  message(FATAL_ERROR "HDF5 >= 1.8.2 or higher has not been found")
endif(HDF5_FOUND)

if(HDF5_IS_PARALLEL)
  set(ALPS_HAVE_HDF5_PARALLEL 1)
  message(
    WARNING
      "parallel(MPI) hdf5 is detected. We will compile but ALPS does not use parallel HDF5. The standard version is preferred."
  )
  if(NOT MPI_FOUND)
    message(
      FATAL_ERROR
        "parallel(MPI) hdf5 needs MPI. Enable MPI or install serial HDF5 libraries."
    )
  endif(NOT MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(HDF5_IS_PARALLEL)

# ##############################################################################
# configuration files
# ##############################################################################

# This should be covered by configure_file, but for some reason it isn't
# probably because the configure file is included at a wrong place so let's set
# it again explicitly
#
if(LAPACK_64_BIT)
  add_definitions(-DBIND_FORTRAN_INTEGER_8)
endif()
if(BIND_FORTRAN_LOWERCASE)
  add_definitions(-DBIND_FORTRAN_LOWERCASE)
endif(BIND_FORTRAN_LOWERCASE)

set(ALPS_XML_PATH "lib")
set(ALPS_XML_DIR "${CMAKE_INSTALL_PREFIX}/${ALPS_XML_PATH}/xml")

# order the include directories so that directories which are in the source or
# build tree always come before directories outside the project.
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

# ##############################################################################
# C & C++ Headers. Apparently required for Mac OS X compilation
# ##############################################################################
include(CheckIncludeFile)
check_include_file(rpc/rpc.h ALPS_HAVE_RPC_XDR_H)
check_include_file(stdarg.h ALPS_HAVE_STDARG_H)
check_include_file(sys/stat.h ALPS_HAVE_SYS_STAT_H)
check_include_file(sys/systeminfo.h ALPS_HAVE_SYS_SYSTEMINFO_H)
check_include_file(sys/types.h ALPS_HAVE_SYS_TYPES_H)
check_include_file(unistd.h ALPS_HAVE_UNISTD_H)
check_include_file(windows.h ALPS_HAVE_WINDOWS_H)

include(CheckIncludeFileCXX)
check_include_file_cxx(valarray ALPS_HAVE_VALARRAY)

# configure files
configure_file(src/alps/config.h.in
               ${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src/alps/config.h)
configure_file(src/alps/version.h.in
               ${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src/alps/version.h)

# ##############################################################################
# include directories
# ##############################################################################
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src)
include_directories(src)

# ##############################################################################
# RPATH setting
# ##############################################################################
if(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
else(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_MACOSX_RPATH 1)
endif(APPLE AND (CMAKE_VERSION VERSION_LESS 2.8.12))

# Disabled tests if (ALPS_BUILD_TESTS AND NOT ALPS_BUILD_SOURCE)
# include(ALPSBuildname) enable_testing() include(CTest) add_subdirectory(test)
# endif (ALPS_BUILD_TESTS AND NOT ALPS_BUILD_SOURCE) add_subdirectory(config)
add_subdirectory(src/alps)
