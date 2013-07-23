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
set (AMBIENT_VERSION_MAJOR 2)
set (AMBIENT_VERSION_MINOR 1)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/ambient" CACHE PATH "Ambient install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(!APPLE)
enable_language(Fortran) # needeed for blas/lapack
endif(!APPLE)

########################################################################
#
# Needed paths
#
########################################################################

set (AMBIENT_DIR ${PROJECT_SOURCE_DIR}/src/ambient)
set (AMBIENT_SOURCES_DIR ${PROJECT_SOURCE_DIR} CACHE STRINGS "Location of the ambient source tree")

########################################################################
#
# Options
#
########################################################################

set (BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Path to the Boost installation (or to the Boost source)")
option(ENABLE_REGRESSION_FUNCTIONAL "Build regression functional tests" OFF)
option(ENABLE_REGRESSION_PERFORMANCE "Build regression performance tests" OFF)
option(ENABLE_REGRESSION_SANDBOX "Build regression sandbox" OFF)

########################################################################
#
# Dependencies
#
########################################################################

if(NOT MPI_CXX_FOUND) # automatic
    find_package(MPI)
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

## if ALPS is present, load its environment
if(ALPS_ROOT_DIR)
    set(ALPS_DIR ${ALPS_ROOT_DIR}/share/alps)
    set(PREVENT_ALPS_COMPILERS ON)
    find_package(ALPS NO_SYSTEM_ENVIRONMENT_PATH)
else()
    set(ALPS_DIR $ENV{ALPS_ROOT} CACHE PATH "Path to the Alps include, need for boost numeric")
endif(ALPS_ROOT_DIR)

## set Boost dependencies
set(Boost_requirements unit_test_framework)
if(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)
  list(APPEND Boost_requirements bindings)
  list(APPEND Boost_requirements system)
endif(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)

## MAQUIS find_boost fall back on ALPS pre-compiled Boost, if ALPS is present
include(${PROJECT_SOURCE_DIR}/config/find_boost.cmake)
find_maquis_boost(COMPONENTS ${Boost_requirements})

message(STATUS "Boost include: ${MAQUIS_Boost_INCLUDE_DIRS}")
message(STATUS "Boost libs: ${MAQUIS_Boost_LIBRARIES}")
message(STATUS "Boost unit_test: ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")

include_directories(${MAQUIS_Boost_INCLUDE_DIRS})
list(APPEND AMBIENT_INCLUDE_DIRS ${MAQUIS_Boost_INCLUDE_DIRS})

link_directories(${MAQUIS_Boost_LIBRARY_DIRS})
list(APPEND AMBIENT_LIBRARY_DIRS ${MAQUIS_Boost_LIBRARY_DIRS})

## Additional libraries / macros for tests
if(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)
  list(APPEND AMBIENT_TESTS_LIBRARIES ${MAQUIS_Boost_LIBRARIES})
    
  macro(compile_with_ambient)
    include_directories(${AMBIENT_INCLUDE_DIRS})
    link_directories(${AMBIENT_LIBRARY_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AMBIENT_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${AMBIENT_EXE_FLAGS}")
  endmacro(compile_with_ambient)
endif(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)

if(NOT BLAS_FOUND)
   find_package(BLAS)
   find_package(LAPACK)
endif(NOT BLAS_FOUND)

if(BLAS_mkl_core_LIBRARY) # because find if does not find the serial
    find_path(BLAS_mkl_INCLUDE_DIR mkl.h ${BLAS_mkl_core_LIBRARY}/../../../include)
    find_library(BLAS_mkl_sequential_LIBRARY mkl_sequential ${BLAS_mkl_INCLUDE_DIR}/../lib/intel64/)
endif()

########################################################################
#
# Compilation
#
########################################################################

include_directories ("${MPI_CXX_INCLUDE_PATH}")
include_directories (src)
add_subdirectory(src/ambient) 

########################################################################
#
# Tests
#
########################################################################

if(ENABLE_REGRESSION_FUNCTIONAL OR ENABLE_REGRESSION_PERFORMANCE)
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

if(ENABLE_REGRESSION_SANDBOX)
  add_subdirectory(regression/sandbox)
endif(ENABLE_REGRESSION_SANDBOX)

########################################################################
#
# Configure files
#
########################################################################

list(APPEND AMBIENT_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include ${MPI_CXX_INCLUDE_PATH})
list(APPEND AMBIENT_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib)
list(APPEND AMBIENT_CONFIG_LIBRARIES ambient ${MPI_CXX_LIBRARIES})
list(APPEND AMBIENT_CONFIG_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
list(APPEND AMBIENT_CONFIG_EXE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
list(APPEND AMBIENT_CONFIG_EXE_FLAGS ${MPI_CXX_LINK_FLAGS})

list(APPEND AMBIENT_BUILD_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src ${MPI_CXX_INCLUDE_PATH})
list(APPEND AMBIENT_BUILD_LIBRARY_DIRS "")
list(APPEND AMBIENT_BUILD_LIBRARIES ambient ${MPI_CXX_LIBRARIES})
list(APPEND AMBIENT_BUILD_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
list(APPEND AMBIENT_BUILD_EXE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
list(APPEND AMBIENT_BUILD_EXE_FLAGS ${MPI_CXX_LINK_FLAGS})

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

