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
message(STATUS "Configuration type: " ${CMAKE_BUILD_TYPE})

set(AMBIENT_SOURCES_DIR "${PROJECT_SOURCE_DIR}/src")

########################################################################
#
# Options
#
########################################################################

set(AMBIENT_THREADING "auto" CACHE STRINGS "cilk, openmp, none or auto.")
set(AMBIENT_MPI "auto" CACHE STRINGS "disable, single, funneled, serialized, multiple or auto.")
set(BLAS_LAPACK_SELECTOR "auto" CACHE STRINGS "auto, manual, mkl_sequential, mkl_parallel, veclib")
set(PLASMA_LIB_DIR "/usr/local/lib" CACHE PATH "Location of PLASMA libraries.")
set(PLASMA_LIBRARIES "-lplasma -lcoreblas -llapacke -lquark" CACHE STRINGS "PLASMA libraries")
set(AMBIENT_CXX_FLAGS "-std=c++11" CACHE STRINGS "C++11 support is required, e.g. -std=c++11 and -gcc-name=XX")
option(AMBIENT_REGRESSION_TESTS "Build regression tests" OFF)

########################################################################
#
# Configuration
#
########################################################################

set(AMBIENT_INCLUDE_DIRS)
set(AMBIENT_LIBRARIES)
set(AMBIENT_LIBRARY_DIRS)
set(AMBIENT_THREADING_FLAGS)
set(AMBIENT_MPI_FLAGS)

if(${AMBIENT_THREADING} MATCHES "cilk")
  set(AMBIENT_THREADING_FLAGS "-DAMBIENT_THREADING=CILK")
elseif(${AMBIENT_THREADING} MATCHES "openmp")
  find_package(OpenMP REQUIRED)
  set(AMBIENT_THREADING_FLAGS "-DAMBIENT_THREADING=OPENMP ${OpenMP_CXX_FLAGS}")
elseif(${AMBIENT_THREADING} MATCHES "none")
  set(AMBIENT_THREADING_FLAGS "-DAMBIENT_THREADING=NONE")
elseif(NOT ${AMBIENT_THREADING} MATCHES "auto")
  message(FATAL_ERROR "AMBIENT_THREADING has to be one of 'cilk', 'openmp', 'none' or 'auto'")
endif()
message(STATUS "Ambient threading: ${AMBIENT_THREADING}")

if(NOT MPI_CXX_FOUND AND NOT ${AMBIENT_MPI} MATCHES "disable")
  find_package(MPI)
endif(NOT MPI_CXX_FOUND AND NOT ${AMBIENT_MPI} MATCHES "disable")
if(MPI_CXX_FOUND AND NOT ${AMBIENT_MPI} MATCHES "disable")
  message(STATUS "MPI include: ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "MPI libs: ${MPI_CXX_LIBRARIES}")
  list(APPEND AMBIENT_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
  list(APPEND AMBIENT_LIBRARIES ${MPI_CXX_LIBRARIES})
else(MPI_CXX_FOUND AND NOT ${AMBIENT_MPI} MATCHES "disable")
  message(STATUS "Warning: MPI is disabled")
  set(AMBIENT_MPI "disable")
endif(MPI_CXX_FOUND AND NOT ${AMBIENT_MPI} MATCHES "disable")

if(${AMBIENT_MPI} MATCHES "disable")
  set(AMBIENT_MPI_FLAGS "-DAMBIENT_MPI=MPI_DISABLE")
elseif(${AMBIENT_MPI} MATCHES "single")
  set(AMBIENT_MPI_FLAGS "-DAMBIENT_MPI=MPI_THREAD_SINGLE")
elseif(${AMBIENT_MPI} MATCHES "funneled")
  set(AMBIENT_MPI_FLAGS "-DAMBIENT_MPI=MPI_THREAD_FUNNELED")
elseif(${AMBIENT_MPI} MATCHES "serialized")
  set(AMBIENT_MPI_FLAGS "-DAMBIENT_MPI=MPI_THREAD_SERIALIZED")
elseif(${AMBIENT_MPI} MATCHES "multiple")
  set(AMBIENT_MPI_FLAGS "-DAMBIENT_MPI=MPI_THREAD_MULTIPLE")
elseif(NOT ${AMBIENT_MPI} MATCHES "auto")
  message(FATAL_ERROR "AMBIENT_MPI has to be one of 'disable', 'single', 'funneled', 'serialized', 'multiple' or 'auto'")
endif()
message(STATUS "Ambient MPI: ${AMBIENT_MPI}")

find_package(MAQUISLapack)
message(STATUS "Math lib dirs: ${MAQUISLapack_LIB_DIRS} ${PLASMA_LIB_DIR}")
message(STATUS "Math libs: ${MAQUISLapack_LIBRARIES} ${PLASMA_LIBRARIES}")
list(APPEND AMBIENT_LIBRARY_DIRS ${MAQUISLapack_LIB_DIRS} ${PLASMA_LIB_DIR})
list(APPEND AMBIENT_LIBRARIES ${MAQUISLapack_LIBRARIES} ${PLASMA_LIBRARIES})

include_directories(${AMBIENT_INCLUDE_DIRS} ${AMBIENT_SOURCES_DIR})
link_directories(${AMBIENT_LIBRARY_DIRS})

########################################################################
#
# Tests
#
########################################################################

if(AMBIENT_REGRESSION_TESTS)
  include(CTest)
  enable_testing() 
  add_subdirectory(tests/functional)
  add_subdirectory(tests/performance)
  if(AMBIENT_THIRD_PARTY_TESTS)
    add_subdirectory(tests/third-party)
  endif(AMBIENT_THIRD_PARTY_TESTS)
else(AMBIENT_REGRESSION_TESTS)
  message("\n\n Ambient is a header-only library (so go ahead and 'make install') \n\n")
endif(AMBIENT_REGRESSION_TESTS)

########################################################################
#
# Configuration dump
#
########################################################################

list(APPEND AMBIENT_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include ${AMBIENT_INCLUDE_DIRS})
list(APPEND AMBIENT_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib ${AMBIENT_LIBRARY_DIRS})
list(APPEND AMBIENT_CONFIG_LIBRARIES ${AMBIENT_LIBRARIES})
list(APPEND AMBIENT_CONFIG_COMPILE_FLAGS ${AMBIENT_CXX_FLAGS} ${AMBIENT_THREADING_FLAGS} ${AMBIENT_MPI_FLAGS})

list(APPEND AMBIENT_BUILD_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src ${AMBIENT_INCLUDE_DIRS})
list(APPEND AMBIENT_BUILD_LIBRARY_DIRS ${AMBIENT_LIBRARY_DIRS})
list(APPEND AMBIENT_BUILD_LIBRARIES ${AMBIENT_LIBRARIES})
list(APPEND AMBIENT_BUILD_COMPILE_FLAGS ${AMBIENT_CXX_FLAGS} ${AMBIENT_THREADING_FLAGS} ${AMBIENT_MPI_FLAGS})

configure_file(AmbientConfig.cmake.in ${PROJECT_BINARY_DIR}/AmbientConfig.cmake)
configure_file(AmbientBuildTree.cmake.in ${PROJECT_BINARY_DIR}/AmbientBuildTree.cmake)

########################################################################
#
# Install
#
########################################################################

install(DIRECTORY src/ambient src/utils
        DESTINATION include
        COMPONENT headers
        PATTERN .svn EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE)

install(DIRECTORY examples
        DESTINATION share
        PATTERN .svn EXCLUDE
        PATTERN CMakeLists.txt EXCLUDE)

install(FILES ${PROJECT_BINARY_DIR}/AmbientConfig.cmake
        DESTINATION share)
