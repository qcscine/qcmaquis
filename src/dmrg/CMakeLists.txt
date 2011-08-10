
cmake_minimum_required(VERSION 2.8)

######################################################################
# Machine Config
######################################################################
set(DEFAULT_BLAS_LAPACK alps)

if(MACHINE_CONFIG)
  if(EXISTS ${MACHINE_CONFIG})
    message(STATUS "Loading config in " ${MACHINE_CONFIG})
    include(${MACHINE_CONFIG})
  else(EXISTS ${MACHINE_CONFIG})
    message(ERROR " Machine config not found!")
  endif(EXISTS ${MACHINE_CONFIG})
endif(MACHINE_CONFIG)

######################################################################
# Project START
######################################################################
project (DMRG) 

# note: you can regenerate operation.pp.h/hpp by invoking smth like the following:
# g++ ambient/core/operation/operation.pp.h.template.h   -P -E -Iboost -I. > operation.pp.h
# g++ ambient/core/operation/operation.pp.hpp.template.h -P -E -Iboost -I. > operation.pp.hpp
# g++ ambient/interface/push.pp.hpp.template.h           -P -E -Iboost -I. > push.pp.hpp

######################################################################
# Libraries
######################################################################
find_package(MPI)

# ALPS
if(NOT ALPS_DIR AND ALPS_ROOT_DIR)
  set(ALPS_DIR ${ALPS_ROOT_DIR}/share/alps)
endif(NOT ALPS_DIR AND ALPS_ROOT_DIR)

set(PREVENT_ALPS_COMPILERS ON)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
#include(${ALPS_USE_FILE})

set(ALPS_LIBRARY_DIRS ${ALPS_LIBRARY_DIRS} CACHE STRINGS "ALPS library dirs")
set(ALPS_INCLUDE_DIRS ${ALPS_INCLUDE_DIRS} CACHE STRINGS "ALPS include dirs")
set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR} CACHE STRINGS "")
set(HDF5_LIBRARIES ${HDF5_LIBRARIES} CACHE STRINGS "")
set(ALPS_LAPACK_DEFINITIONS ${LAPACK_DEFINITIONS} CACHE STRINGS "")
set(ALPS_LAPACK_LINKER_FLAGS " ${LAPACK_LINKER_FLAGS}" CACHE STRINGS "")
set(ALPS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRINGS "")
set(ALPS_LAPACK_LIBRARY ${LAPACK_LIBRARY} CACHE STRINGS "")
set(ALPS_BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE STRINGS "")
set(ALPS_BLAS_LIBRARY ${BLAS_LIBRARY} CACHE STRINGS "")

mark_as_advanced(ALPS_LIBRARY_DIRS ALPS_INCLUDE_DIRS HDF5_INCLUDE_DIR HDF5_LIBRARIES
                ALPS_LAPACK_DEFINITIONS ALPS_LAPACK_LINKER_FLAGS ALPS_LAPACK_LIBRARIES ALPS_LAPACK_LIBRARY
                ALPS_BLAS_LIBRARIES ALPS_BLAS_LIBRARY)

# OpenMP
find_package(OpenMP)

# MPI
if(NOT MPI_FOUND)
  find_package(MPI)
endif(NOT MPI_FOUND)

######################################################################
# Version information
######################################################################
set(DMRG_YEAR 2011)
set(DMRG_VERSION_MAJOR 0)
set(DMRG_VERSION_MINOR 1)
set(DMRG_VERSION_BUILD "")

if(NOT DMRG_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
  FIND_PACKAGE(Subversion) 
  IF(Subversion_FOUND) 
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} DMRG) 
    set(DMRG_VERSION_BUILD "r${DMRG_WC_REVISION}")
  ENDIF(Subversion_FOUND) 
endif(NOT DMRG_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

if(DMRG_VERSION_BUILD)
  set(DMRG_VERSION "${DMRG_VERSION_MAJOR}.${DMRG_VERSION_MINOR}-${DMRG_VERSION_BUILD}")
else(DMRG_VERSION_BUILD)
  set(DMRG_VERSION "${DMRG_VERSION_MAJOR}.${DMRG_VERSION_MINOR}")
endif(DMRG_VERSION_BUILD)
set(DMRG_VERSION_STRING "DMRG version ${DMRG_VERSION}")
MESSAGE(STATUS "DMRG version: ${DMRG_VERSION}") 

######################################################################
# CPU Freq.
######################################################################

# TODO: find a way to get this information during config
# Mac: sysctl -a hw
# Unix: cat /proc/cpuinfo
set(CPU_FREQ 2.8e9 CACHE STRING "CPU Frequency")
add_definitions(-DCPU_FREQ=${CPU_FREQ})


######################################################################
# What to build
######################################################################

#set(CMAKE_C_FLAGS "")
#set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
#set(CMAKE_C_FLAGS_DEBUG "-g")

#set(CMAKE_CXX_FLAGS "")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")

######################################################################
# User options
######################################################################

option(ENABLE_OMP "Enable OpenMP" OFF)
option(ENABLE_AMBIENT "Enable ambient" OFF)
option(USE_MTM "Use multi threads matrix class" OFF)
option(USE_GPU "Use GPU for large matrix gemm (to be used together with USE_MTM)" OFF)
set(EXTRA_DEFINITIONS "-g -O3 -DNDEBUG" CACHE STRING "Advanced definitions")

option(BUILD_APP "Build DMRG Apps" ON)
option(BUILD_BELA_BRANCH "Build Bela's branch DMRG Apps" OFF)
set(BELA_BRANCH ${PROJECT_SOURCE_DIR}/../../branches/bela/src CACHE FILEPATH "Path to Bela's branch")
option(BUILD_REGRESSION "Build regression tests" OFF)
option(BUILD_MULTIGRID "Build multigrid code" OFF)

set(BLAS_LAPACK_SELECTOR ${DEFAULT_BLAS_LAPACK} CACHE STRING "Choose the Blas/Lapack version to use, possible values 'alps', 'veclib', 'manual'")
set(BLAS_LAPACK_MANUAL_LIBS ${BLAS_LAPACK_MANUAL_LIBS} CACHE STRINGS "")
set(BLAS_LAPACK_MANUAL_LIBS_DIR ${BLAS_LAPACK_MANUAL_LIBS_DIR} CACHE STRINGS "")
set(BLAS_LAPACK_MANUAL_DEFINITIONS ${BLAS_LAPACK_MANUAL_DEFINITIONS} CACHE STRINGS "")
set(BLAS_LAPACK_MANUAL_INCLUDES ${BLAS_LAPACK_MANUAL_INCLUDES} CACHE STRINGS "")

######################################################################
# Parsing user options
######################################################################

set(DMRG_LIBRARIES alps boost ${HDF5_LIBRARIES})

# BLAS/LAPACK selector
if(${BLAS_LAPACK_SELECTOR} MATCHES "alps")
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${ALPS_LAPACK_LINKER_FLAGS})
  set(BLAS_LAPACK_LIBRARIES ${ALPS_LAPACK_LIBRARIES} ${ALPS_LAPACK_LIBRARY}
      ${ALPS_BLAS_LIBRARIES} ${ALPS_BLAS_LIBRARY})
  set(DMRG_LIBRARIES ${DMRG_LIBRARIES} ${BLAS_LAPACK_LIBRARIES})
  add_definitions(${ALPS_LAPACK_DEFINITIONS})
endif(${BLAS_LAPACK_SELECTOR} MATCHES "alps")
if(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " -framework vecLib")
    set(BLAS_LAPACK_LIBRARIES /usr/lib/libpthread.dylib)
    set(DMRG_LIBRARIES ${DMRG_LIBRARIES} ${BLAS_LAPACK_LIBRARIES})
  else(APPLE)
    meassage("Error: VecLib available only on Mac.")
  endif(APPLE)
endif(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
if(${BLAS_LAPACK_SELECTOR} MATCHES "manual")
  set(BLAS_LAPACK_LIBRARIES ${BLAS_LAPACK_MANUAL_LIBS})
  set(DMRG_LIBRARIES ${DMRG_LIBRARIES} ${BLAS_LAPACK_LIBRARIES})
  add_definitions(${BLAS_LAPACK_MANUAL_DEFINITIONS})
  include_directories(${BLAS_LAPACK_MANUAL_INCLUDES})
  link_directories(${BLAS_LAPACK_MANUAL_LIBS_DIR})
endif(${BLAS_LAPACK_SELECTOR} MATCHES "manual")

# Extra definitions
add_definitions(${EXTRA_DEFINITIONS})
if(USE_MTM)
  add_definitions(-DUSE_MTM)
endif(USE_MTM)

# Include/Link Directories
include_directories(${PROJECT_SOURCE_DIR} ${HDF5_INCLUDE_DIR} ${ALPS_INCLUDE_DIRS})
link_directories(${ALPS_LIBRARY_DIRS})

# OpenMPI
if(ENABLE_OMP)
  if(OPENMP_FOUND)
    macro(enable_omp_if_found)
      #add_definitions(${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
      set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
    endmacro(enable_omp_if_found)
  else(OPENMP_FOUND)
    message("OpenMP not found")
  endif(OPENMP_FOUND)
else(ENABLE_OMP)
  macro(enable_omp_if_found)
  endmacro(enable_omp_if_found)
endif(ENABLE_OMP)

# Ambient & MPI
if(MPI_FOUND)
  macro(enable_mpi)
    add_definitions(${MPI_COMPILE_FLAGS})
    #set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS})
    #set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
  endmacro(enable_mpi)
else(MPI_FOUND)
  macro(enable_mpi)
    message("MPI not found!")
  endmacro(enable_mpi)
endif(MPI_FOUND)
if(ENABLE_AMBIENT)
  if (MPI_FOUND)
    macro(enable_ambient)
      enable_mpi()
      add_definitions(-DMPI_PARALLEL)
      set(AMBIENT_LIBRARIES ambient boost ${MPI_LIBRARIES})
      set(AMBIENT_SCALAPACK_LIBRARIES -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64 ${AMBIENT_LIBRARIES})
    endmacro(enable_ambient)
  else(MPI_FOUND)
    message("MPI not found!")
  endif(MPI_FOUND)
else(ENABLE_AMBIENT)
  macro(enable_ambient)
    message("Ambient not enabled!")
  endmacro(enable_ambient)
endif(ENABLE_AMBIENT)

# GPU / CUDA
# TODO

######################################################################
# Targets
######################################################################

if(ENABLE_AMBIENT)
  add_subdirectory(ambient)
endif(ENABLE_AMBIENT)
if(BUILD_APP)
  add_subdirectory(app)
else(BUILD_APP)
  add_subdirectory(app ECLUDE_FROM_ALL)
endif(BUILD_APP)
if(BUILD_BELA_BRANCH)
  if (NOT EXISTS ${PROJECT_SOURCE_DIR}/belas_branch)
    execute_process(COMMAND ln -s ${BELA_BRANCH} ${PROJECT_SOURCE_DIR}/belas_branch)
  endif(NOT EXISTS ${PROJECT_SOURCE_DIR}/belas_branch)
  add_subdirectory(belas_branch)
endif(BUILD_BELA_BRANCH)

if(BUILD_REGRESSION)
  add_subdirectory(regression/dmrg)
  if(ENABLE_AMBIENT)
    add_subdirectory(regression/p_dense_matrix)
    add_subdirectory(regression/p_memory_management)
    add_subdirectory(regression/p_diagonal_matrix)
    add_subdirectory(regression/p_block_matrix)
    add_subdirectory(regression/p_dense_matrix_benchmark)
    add_subdirectory(regression/pdmrg)
  endif(ENABLE_AMBIENT)
else(BUILD_REGRESSION)
  add_subdirectory(regression/dmrg EXCLUDE_FROM_ALL)
  if(ENABLE_AMBIENT)
    add_subdirectory(regression/p_dense_matrix EXCLUDE_FROM_ALL)
    add_subdirectory(regression/p_memory_management EXCLUDE_FROM_ALL)
    add_subdirectory(regression/p_diagonal_matrix EXCLUDE_FROM_ALL)
    add_subdirectory(regression/p_block_matrix EXCLUDE_FROM_ALL)
    add_subdirectory(regression/p_dense_matrix_benchmark EXCLUDE_FROM_ALL)
    add_subdirectory(regression/pdmrg EXCLUDE_FROM_ALL)
  endif(ENABLE_AMBIENT)
endif(BUILD_REGRESSION)
