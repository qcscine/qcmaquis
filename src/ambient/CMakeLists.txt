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
if(NOT MPI_CXX_FOUND) # automatic
    find_package(MPI)
else() # by hand, I only set up the new flags of FindMPI to respect the new norm
    set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} CACHE STRINGS "MPI CXX include path")
    set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} CACHE STRINGS "MPI CXX linking flags" )
    set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} CACHE STRINGS "MPI CXX libraries to link against")
    set(MPIEXEC ${MPIEXEC} CACHE STRINGS "Execute MPI applications")
    set(MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG} CACHE STRINGS "Flag used by MPI to specify the number of processes for MPIEXEC; the next option will be the number of processes.")
    set(MPIEXEC_PREFLAGS ${MPIEXEC_PREFLAGS} CACHE STRINGS "These flags will be directly before the executable that is being run by MPIEXEC.")
    set(MPIEXEC_POSTFLAGS ${MPIEXEC_POSTFLAGS} CACHE STRINGS "These flags will come after all flags given to MPIEXEC.")
    set(MPIEXEC_MAX_NUMPROCS ${MPIEXEC_MAX_NUMPROCS} "2" CACHE STRING "Maximum number of processors available to run MPI applications.")
endif(NOT MPI_CXX_FOUND)

########################################################################
#
# Options
#
########################################################################


set (BOOST_ROOT $ENV{BOOST_ROOT} CACHE PATH "Path to the Boost installation (or to the Boost source)")
option (AMBIENT_DEBUG_COMPILE "Create a debug compile of AMBIENT" OFF)
option (AMBIENT_MAIN_COMPILE "Compile the sandbox " OFF)
option (AMBIENT_TESTS "Build and execute the AMBIENT regression tests" OFF)
option (BLAS_MULTITHREADED "Build and execute BLAS-Lapack multithreaded" OFF)
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
include_directories ("${MPI_CXX_INCLUDE_PATH}")

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

if(NOT BLAS_FOUND)
   find_package(BLAS)
   find_package(LAPACK)
endif(NOT BLAS_FOUND)

if(BLAS_mkl_core_LIBRARY) #because find if does not find the serial
    find_path(BLAS_mkl_INCLUDE_DIR mkl.h ${BLAS_mkl_core_LIBRARY}/../../../include)
    find_library(BLAS_mkl_sequential_LIBRARY mkl_sequential ${BLAS_mkl_INCLUDE_DIR}/../lib/intel64/)
endif()


########################################################################
#
# Main Compile
#
########################################################################

if(AMBIENT_MAIN_COMPILE)
    add_subdirectory(sandbox)
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

