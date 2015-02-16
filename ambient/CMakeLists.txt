cmake_minimum_required (VERSION 2.8.5)

########################################################################
#
# Project information
#
########################################################################

if(MACHINE_CONFIG)
  if(EXISTS ${MACHINE_CONFIG})
    message(STATUS "Loading config in " ${MACHINE_CONFIG})
    include(${MACHINE_CONFIG})
  else(EXISTS ${MACHINE_CONFIG})
    message(ERROR " Machine config not found!")
  endif(EXISTS ${MACHINE_CONFIG})
endif(MACHINE_CONFIG)

project (ambient)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/ambient" CACHE PATH "Ambient install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

set(AMBIENT_SOURCES_DIR "${PROJECT_SOURCE_DIR}/src")

########################################################################
#
# Options
#
########################################################################

set(AMBIENT_THREADING "auto" CACHE STRINGS "Ambient threading backend: cilk, openmp, serial or auto.")
set(BLAS_LAPACK_SELECTOR "auto" CACHE STRINGS "auto, manual, mkl_sequential, mkl_parallel, veclib")
set(PLASMA_LIB_DIR "/usr/local/lib" CACHE PATH "Location of PLASMA libraries.")
set(PLASMA_LIBRARIES "-lplasma -lcoreblas -llapacke -lquark" CACHE STRINGS "PLASMA libraries")
set(AMBIENT_CXX_FLAGS "-std=c++11" CACHE STRINGS "C++11 support is required, e.g. -std=c++11 and -gcc-name=XX")

option(AMBIENT_DISABLE_MPI "Disable MPI support (shared memory only)" OFF)
option(FUNCTIONAL_REGRESSION_TESTS  "Build functional regression tests" OFF)
option(PERFORMANCE_REGRESSION_TESTS "Build performance regression tests" OFF)
option(THIRD_PARTY_REGRESSION_TESTS "Build third-party regression tests" OFF)

########################################################################
#
# Dependencies
#
########################################################################

set(AMBIENT_INCLUDE_DIRS)
set(AMBIENT_LIBRARIES)
set(AMBIENT_LIBRARY_DIRS)

if(${AMBIENT_THREADING} MATCHES "cilk")
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_CILK")
elseif(${AMBIENT_THREADING} MATCHES "openmp")
  find_package(OpenMP REQUIRED)
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_OMP ${OpenMP_CXX_FLAGS}")
elseif(${AMBIENT_THREADING} MATCHES "serial")
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_SERIAL")
elseif(NOT ${AMBIENT_THREADING} MATCHES "auto")
  message(FATAL_ERROR "Variable AMBIENT_THREADING has to be one of 'cilk', 'openmp', 'serial' or 'auto'")
endif()
message(STATUS "Ambient threading: ${AMBIENT_THREADING}")

if(AMBIENT_DISABLE_MPI)
  message(STATUS "MPI is disabled")
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_DISABLE_MPI")
else(AMBIENT_DISABLE_MPI)
  if(NOT MPI_CXX_FOUND)
    find_package(MPI REQUIRED)
  else()
    set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} CACHE STRINGS "MPI CXX include path")
    set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} CACHE STRINGS "MPI CXX linking flags" )
    set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} CACHE STRINGS "MPI CXX libraries to link against")
  endif(NOT MPI_CXX_FOUND)

  message(STATUS "MPI include: ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "MPI libs: ${MPI_CXX_LIBRARIES}")
  list(APPEND AMBIENT_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
  list(APPEND AMBIENT_LIBRARIES ${MPI_CXX_LIBRARIES})
endif(AMBIENT_DISABLE_MPI)

find_package(MAQUISLapack)
message(STATUS "Math lib dirs: ${MAQUISLapack_LIB_DIRS} ${PLASMA_LIB_DIR}")
message(STATUS "Math libs:     ${MAQUISLapack_LIBRARIES} ${PLASMA_LIBRARIES}")
list(APPEND AMBIENT_LIBRARY_DIRS ${MAQUISLapack_LIB_DIRS} ${PLASMA_LIB_DIR})
list(APPEND AMBIENT_LIBRARIES ${MAQUISLapack_LIBRARIES} ${PLASMA_LIBRARIES})

########################################################################
#
# Compilation
#
########################################################################

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMBIENT_CXX_FLAGS}")
include_directories(${AMBIENT_INCLUDE_DIRS} ${AMBIENT_SOURCES_DIR})
link_directories(${AMBIENT_LIBRARY_DIRS})
add_subdirectory(src/ambient)

########################################################################
#
# Tests
#
########################################################################

if(THIRD_PARTY_REGRESSION_TESTS)

  if(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})
    set(ALPS_DIR $ENV{ALPS_ROOT}/share/alps)
    set(PREVENT_ALPS_COMPILERS ON)
    find_package(ALPS NO_SYSTEM_ENVIRONMENT_PATH)
  endif(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})

  include(${PROJECT_SOURCE_DIR}/config/find_boost.cmake)
  find_maquis_boost(COMPONENTS)
  message(STATUS "Boost include: ${MAQUIS_Boost_INCLUDE_DIRS}")
  message(STATUS "Boost libs: ${MAQUIS_Boost_LIBRARIES}")

  list(APPEND AMBIENT_INCLUDE_DIRS ${MAQUIS_Boost_INCLUDE_DIRS})
  list(APPEND AMBIENT_THIRD_PARTY_LIBRARIES ${AMBIENT_LIBRARIES} -lboost_system)
  
endif(THIRD_PARTY_REGRESSION_TESTS)
if(FUNCTIONAL_REGRESSION_TESTS OR PERFORMANCE_REGRESSION_TESTS OR THIRD_PARTY_REGRESSION_TESTS)

  macro(compile_with_ambient)
    include_directories(${AMBIENT_INCLUDE_DIRS})
    link_directories(${AMBIENT_LIBRARY_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMBIENT_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${AMBIENT_EXE_FLAGS}")
  endmacro(compile_with_ambient)

  include(CTest)
  enable_testing() 

  if(FUNCTIONAL_REGRESSION_TESTS)
    add_subdirectory(tests/functional)
  endif(FUNCTIONAL_REGRESSION_TESTS)

  if(PERFORMANCE_REGRESSION_TESTS)
    add_subdirectory(tests/performance)
  endif(PERFORMANCE_REGRESSION_TESTS)

  if(THIRD_PARTY_REGRESSION_TESTS)
    add_subdirectory(tests/third-party)
  endif(THIRD_PARTY_REGRESSION_TESTS)

endif(FUNCTIONAL_REGRESSION_TESTS OR PERFORMANCE_REGRESSION_TESTS OR THIRD_PARTY_REGRESSION_TESTS)

########################################################################
#
# Configuration dump
#
########################################################################

list(APPEND AMBIENT_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include ${AMBIENT_INCLUDE_DIRS})
list(APPEND AMBIENT_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib ${AMBIENT_LIBRARY_DIRS})
list(APPEND AMBIENT_CONFIG_LIBRARIES ambient ${AMBIENT_LIBRARIES})
list(APPEND AMBIENT_CONFIG_COMPILE_FLAGS ${AMBIENT_CXX_FLAGS})

list(APPEND AMBIENT_BUILD_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src ${AMBIENT_INCLUDE_DIRS})
list(APPEND AMBIENT_BUILD_LIBRARY_DIRS ${AMBIENT_LIBRARY_DIRS})
list(APPEND AMBIENT_BUILD_LIBRARIES ambient ${AMBIENT_LIBRARIES})
list(APPEND AMBIENT_BUILD_COMPILE_FLAGS ${AMBIENT_CXX_FLAGS})

configure_file(AmbientConfig.cmake.in ${PROJECT_BINARY_DIR}/AmbientConfig.cmake)
configure_file(AmbientBuildTree.cmake.in ${PROJECT_BINARY_DIR}/AmbientBuildTree.cmake)

########################################################################
#
# Install
#
########################################################################

install(FILES ${PROJECT_BINARY_DIR}/AmbientConfig.cmake DESTINATION share)

install(DIRECTORY src/ambient src/utils
        DESTINATION include
        COMPONENT headers
        PATTERN .svn EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE
       )
