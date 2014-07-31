cmake_minimum_required (VERSION 2.8.5)

########################################################################
#
# Project and version information
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
set (AMBIENT_YEAR 2014)
set (AMBIENT_VERSION_MAJOR 2)
set (AMBIENT_VERSION_MINOR 1)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/ambient" CACHE PATH "Ambient install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(AMBIENT_SOURCES_DIR "${PROJECT_SOURCE_DIR}/src")

if(NOT APPLE)
  enable_language(Fortran) # needeed for blas/lapack
endif(NOT APPLE)


########################################################################
#
# Options
#
########################################################################

set(AMBIENT_THREADING "auto" CACHE STRINGS "Ambient threading backend: cilk, openmp, serial or auto.")
set(BLAS_LAPACK_SELECTOR "auto" CACHE STRINGS "auto, manual, mkl_sequential, mkl_parallel, veclib, alps")
set(PLASMA_LIB_DIR "/usr/local/lib" CACHE PATH "Location of PLASMA libraries.")
set(PLASMA_LIBRARIES "-lplasma -lcoreblas -llapacke -lquark" CACHE STRINGS "PLASMA libraries")
set(AMBIENT_CPP11 "-std=c++11" CACHE STRINGS "Flags needed for C++11, e.g. -std=c++11 and -gcc-name=XX. Leave empty to disable C++11.")

set (BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Path to the Boost installation (or to the Boost source)")
option(ENABLE_REGRESSION_FUNCTIONAL "Build regression functional tests" OFF)
option(ENABLE_REGRESSION_PERFORMANCE "Build regression performance tests" OFF)


########################################################################
#
# Find dependencies
#
########################################################################

set(AMBIENT_INCLUDE_DIRS)
set(AMBIENT_LIBRARIES)
set(AMBIENT_LIBRARY_DIRS)
set(AMBIENT_CXX_FLAGS)

# MPI
#####
if(NOT MPI_CXX_FOUND) # automatic
  find_package(MPI REQUIRED)
else() # by hand, I only set up the new flags of FindMPI to respect the new norm
  set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} CACHE STRINGS "MPI CXX include path")
  set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} CACHE STRINGS "MPI CXX linking flags" )
  set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} CACHE STRINGS "MPI CXX libraries to link against")
  set(MPIEXEC ${MPIEXEC} CACHE STRINGS "Execute MPI applications")
  set(MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG} CACHE STRINGS "Flag used by MPI to specify the number of processes for MPIEXEC.")
  set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} CACHE STRINGS "These flags will be directly before the executable that is being run by MPIEXEC.")
  set(MPIEXEC_POSTFLAGS ${MPIEXEC_POSTFLAGS} CACHE STRINGS "These flags will come after all flags given to MPIEXEC.")
  set(MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS} "2" CACHE STRING "Maximum number of processors available to run MPI applications.")
endif(NOT MPI_CXX_FOUND)

message(STATUS "MPI include: ${MPI_CXX_INCLUDE_PATH}")
message(STATUS "MPI libs: ${MPI_CXX_LIBRARIES}")
list(APPEND AMBIENT_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
list(APPEND AMBIENT_LIBRARIES ${MPI_CXX_LIBRARIES})

# Threading
###########
if(${AMBIENT_THREADING} MATCHES "cilk")
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_CILK")
elseif(${AMBIENT_THREADING} MATCHES "openmp")
  find_package(OpenMP REQUIRED)
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_OMP ${OpenMP_CXX_FLAGS}")
elseif(${AMBIENT_THREADING} MATCHES "serial")
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} -DAMBIENT_SERIAL")
elseif(NOT ${AMBIENT_THREADING} MATCHES "auto")
  message(FATAL_ERROR "Variable AMBIENT_THREADING has to be one of 'cilk', 'openmp', 'serial' or 'auto'.")
endif()
message(STATUS "Ambient threading: ${AMBIENT_THREADING}")

# CPP11
#######
if(NOT "${AMBIENT_CPP11}" STREQUAL "")
  set(AMBIENT_CXX_FLAGS "${AMBIENT_CXX_FLAGS} ${AMBIENT_CPP11}")
endif(NOT "${AMBIENT_CPP11}" STREQUAL "")

# ALPS & Boost
##############
if(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})
  set(ALPS_DIR $ENV{ALPS_ROOT}/share/alps)
  set(PREVENT_ALPS_COMPILERS ON)
  find_package(ALPS NO_SYSTEM_ENVIRONMENT_PATH)
endif(NOT ALPS_DIR AND DEFINED ENV{ALPS_ROOT})

## set Boost dependencies
set(Boost_requirements)
if(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)
  list(APPEND Boost_requirements bindings)
  list(APPEND Boost_requirements system)
endif(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)

## MAQUIS find_boost fall back on ALPS pre-compiled Boost, if ALPS is present
include(${PROJECT_SOURCE_DIR}/config/find_boost.cmake)
find_maquis_boost(COMPONENTS ${Boost_requirements})

message(STATUS "Boost include: ${MAQUIS_Boost_INCLUDE_DIRS}")
message(STATUS "Boost libs: ${MAQUIS_Boost_LIBRARIES}")
list(APPEND AMBIENT_INCLUDE_DIRS ${MAQUIS_Boost_INCLUDE_DIRS})
list(APPEND AMBIENT_LIBRARY_DIRS ${MAQUIS_Boost_LIBRARY_DIRS})

# BLAS & LAPACK
###############
find_package(MAQUISLapack)
message(STATUS "Lapack include:  ${MAQUISLapack_INCLUDE_DIRS}")
message(STATUS "Lapack lib dirs: ${MAQUISLapack_LIB_DIRS} ${PLASMA_LIB_DIR}")
message(STATUS "Lapack libs:     ${MAQUISLapack_LIBRARIES} ${PLASMA_LIBRARIES}")
list(APPEND AMBIENT_INCLUDE_DIRS ${MAQUISLapack_INCLUDE_DIRS})
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

if(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)
  list(APPEND AMBIENT_TESTS_LIBRARIES ambient ${AMBIENT_LIBRARIES} ${MAQUIS_Boost_LIBRARIES})
  
  macro(compile_with_ambient)
    include_directories(${AMBIENT_INCLUDE_DIRS})
    link_directories(${AMBIENT_LIBRARY_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMBIENT_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${AMBIENT_EXE_FLAGS}")
  endmacro(compile_with_ambient)

  include(CTest)
  enable_testing() 
  include_directories(${PROJECT_SOURCE_DIR}/regression)
endif(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)

if(ENABLE_REGRESSION_FUNCTIONAL)
  add_subdirectory(regression/functional)
endif(ENABLE_REGRESSION_FUNCTIONAL)

if(ENABLE_REGRESSION_PERFORMANCE)
  add_subdirectory(regression/performance)
endif(ENABLE_REGRESSION_PERFORMANCE)

########################################################################
#
# Configure files
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

#install(FILES ${PROJECT_BINARY_DIR}/ambient/libambient.so DESTINATION lib)

install(FILES ${PROJECT_BINARY_DIR}/AmbientConfig.cmake DESTINATION share)
install(EXPORT AmbientTargets DESTINATION share COMPONENT dev)

set(AMBIENT_INSTALL_HEADERS src/ambient/ambient.h
                            src/ambient/ambient.hpp
                            src/ambient/auxiliary.h
                            src/ambient/model.h
                            src/ambient/groups/group.h
                            src/ambient/groups/multirank.h
                            src/ambient/packets/types.h
                            src/ambient/packets/packet.h
                            src/ambient/packets/packet_t.h
                            src/ambient/core/auxiliary.h
                            src/ambient/core/layout.h
                            src/ambient/core/operation/operation.h
                            src/ambient/core/operation/operation.pp.hpp
                            src/ambient/core/memblock.h
                            src/ambient/core/p_object.h
                            src/ambient/core/scope_context.h
                            src/ambient/core/select.h
                            src/ambient/core/workgroup_context.h
                            src/ambient/interface/core.hpp
                            src/ambient/interface/c_kernels.hpp
                            src/ambient/interface/interface.hpp
                            src/ambient/interface/l_kernels.hpp
                            src/ambient/interface/model.hpp
                            src/ambient/interface/push.pp.hpp
                            src/ambient/groups/packet_manager.h)

foreach(F ${AMBIENT_INSTALL_HEADERS})
  list(APPEND INSTALL_REGEX PATTERN ${F})
endforeach(F)

install(DIRECTORY src/ambient src/utils
        DESTINATION include
        COMPONENT headers
        FILES_MATCHING ${INSTALL_REGEX}
        PATTERN .svn EXCLUDE
        )

