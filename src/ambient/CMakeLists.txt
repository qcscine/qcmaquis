cmake_minimum_required (VERSION 2.8)

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
set (AMBIENT_VERSION_MAJOR 0)
set (AMBIENT_VERSION_MINOR 1)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/maquis" CACHE PATH "MAQUIS install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if(!APPLE)
enable_language(Fortran) #needeed for blas/lapack
endif(!APPLE)
########################################################################
#
# several needed  paths 
#
########################################################################

set (AMBIENT_DIR ${PROJECT_SOURCE_DIR}/src/ambient)
set (P_DENSE_MATRIX_DIR ${PROJECT_SOURCE_DIR}/../types/src)

########################################################################
#
# Compilers mpi
#
########################################################################
if(NOT MPI_FOUND) # automatic
    find_package(MPI)
#    include(CMakeForceCompiler)
#    CMAKE_FORCE_C_COMPILER(mpiicc ${CMAKE_C_COMPILER_ID})
#    CMAKE_FORCE_CXX_COMPILER(mpiicpc  ${CMAKE_CXX_COMPILER_ID})
endif(NOT MPI_FOUND)

########################################################################
#
# Options
#
########################################################################

if(MPI_FOUND) # we used a config file
set(MPI_COMPILE_FLAGS ${MPI_COMPILE_FLAGS} CACHE STRINGS "")
set(MPI_INCLUDE_PATH ${MPI_INCLUDE_PATH} CACHE STRINGS "" )
set(MPI_LINK_FLAGS ${MPI_LINK_FLAGS} CACHE STRINGS "" )
set(MPI_LIBRARY ${MPI_LIBRARY} CACHE STRINGS "")
set(MPI_EXTRA_LIBRARY ${MPI_EXTRA_LIBRARY} CACHE STRINGS "")
set(MPI_CXX_INCLUDE_PATH ${MPI_INCLUDE_PATH})
set(MPI_CXX_LIBRARIES ${MPI_LIBRARY})
set(MPIEXEC ${MPIEXEC})
set(MPIEXEC_NUMPROC_FLAG ${PMPIEXEC_NUMPROC_FLAG} CACHE STRINGS "")
set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} CACHE STRINGS "")
set(MPIEXEC_POSTFLAGS ${MPIEXEC_POSTFLAGS} CACHE STRINGS "")
endif(MPI_FOUND)

set (BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Path to the Boost installation (or to the Boost source)")
option (AMBIENT_DEBUG_COMPILE "Create a debug compile of AMBIENT" OFF)
option (AMBIENT_MAIN_COMPILE "Compile the sandbox " OFF)
option (AMBIENT_TESTS "Build and execute the AMBIENT regression tests" OFF)
option (BLAS_MULTITHREADED "Build and execute BLAS-Lapack multithreaded" ON)
######################################################################
# CPU Freq.
######################################################################

# TODO: find a way to get this information during config
# C - Dangerous the CPU frequence is sometimes between front and compute nodes
# Mac: sysctl -a hw
# Unix: cat /proc/cpuinfo
set(CPU_FREQ 2.8e9 CACHE STRING "CPU Frequency")
add_definitions(-DCPU_FREQ=${CPU_FREQ})

########################################################################
#
# Find dependencies
#
########################################################################

#set(BOOST_ROOT ${BOOST_DIR})
find_package (Boost COMPONENTS unit_test_framework)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost libraries not found. Please specify location using the BOOST_DIR variable")
endif()

########################################################################
#
# Ambient compile 
#
########################################################################
#include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${MPI_INCLUDE_PATH}")

include_directories (src)
add_subdirectory(src/ambient) 

########################################################################
#
# Library BLAS/LAPACK
#
# specify your blas : MKL = Intel10_64lp
#              acml       = ACML
#
########################################################################

#set(BLA_VENDOR veclib CACHE STRING "BLAS VENDOR")

if(NOT BLAS_FOUND)
   find_package(BLAS)
   find_package(LAPACK)

#   if(${BLA_VENDOR} STREQUAL "Intel10_64lp")
#       find_path(MKL_INCLUDE_DIR mkl.h $ENV{MKLROOT}/include)
#       include_directories("${MKL_INCLUDE_DIR}")

#Find SCaLAPACK and BLACS FOR OPENMPI
#       find_library(SCALAPACK_LIBRARY mkl_scalapack_lp64  $ENV{MKLROOT}/lib/em64t)
#       find_library(BLACS_LIBRARY mkl_blacs_openmpi_lp64 $ENV{MKLROOT}/lib/em64t)
#       find_library(MKL_SEQUENTIAL mkl_sequential $ENV{MKLROOT}/lib/em64t)

#   endif()
endif(NOT BLAS_FOUND)

########################################################################
#
# Main Compile
#
########################################################################

if(AMBIENT_MAIN_COMPILE)
    add_subdirectory(sandbox)
#    if(${BLA_VENDOR} STREQUAL "Intel10_64lp")
#        target_link_libraries(main.out ambient ${SCALAPACK_LIBRARY} ${BLAS_mkl_intel_lp64_LIBRARY} ${BLAS_mkl_core_LIBRARY} ${BLACS_LIBRARY} ${MKL_SEQUENTIAL} pthread)
#    endif()
endif()

########################################################################
#
# Tests
#
########################################################################

if(AMBIENT_TESTS)
    enable_testing()
    add_subdirectory (regression/)
endif(AMBIENT_TESTS)


########################################################################
#
# Configure files
#
########################################################################
list(APPEND AMBIENT_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include ${MPI_CXX_INCLUDE_PATH})
list(APPEND AMBIENT_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib)
list(APPEND AMBIENT_CONFIG_LIBRARIES ambient ${MPI_CXX_LIBRARIES})
list(APPEND AMBIENT_CONFIG_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
list(APPEND AMBIENT_CONFIG_EXE_FLAGS -DMPI_PARALLEL ${MPI_CXX_COMPILE_FLAGS})
# list(APPEND AMBIENT_CONFIG_EXE_FLAGS ${MPI_CXX_LINK_FLAGS})

list(APPEND AMBIENT_BUILD_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src ${MPI_CXX_INCLUDE_PATH})
list(APPEND AMBIENT_BUILD_LIBRARY_DIRS "")
list(APPEND AMBIENT_BUILD_LIBRARIES ambient ${MPI_CXX_LIBRARIES})
list(APPEND AMBIENT_BUILD_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS})
list(APPEND AMBIENT_BUILD_EXE_FLAGS -DMPI_PARALLEL ${MPI_CXX_COMPILE_FLAGS})
# list(APPEND AMBIENT_BUILD_EXE_FLAGS ${MPI_CXX_LINK_FLAGS})

configure_file(MaquisAmbientConfig.cmake.in ${PROJECT_BINARY_DIR}/MaquisAmbientConfig.cmake)
configure_file(MaquisAmbientBuildTree.cmake.in ${PROJECT_BINARY_DIR}/MaquisAmbientBuildTree.cmake)


########################################################################
#
# Install
#
########################################################################

#install(FILES ${PROJECT_BINARY_DIR}/ambient/libambient.so DESTINATION lib)

install(FILES ${PROJECT_BINARY_DIR}/MaquisAmbientConfig.cmake DESTINATION share)

install(EXPORT MaquisAmbientTargets DESTINATION share COMPONENT dev)

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

