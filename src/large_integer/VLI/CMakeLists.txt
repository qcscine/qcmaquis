cmake_minimum_required (VERSION 2.6)


########################################################################
#
# DEBUG move false to true
#
########################################################################
set(CMAKE_VERBOSE_MAKEFILE false)

########################################################################
#
# Project and version information
#
########################################################################

project (vli)
set (VLI_VERSION_MAJOR 0)
set (VLI_VERSION_MINOR 1)



########################################################################
#
# Options
#
########################################################################

set (BOOST_DIR $ENV{BOOST_DIR} CACHE PATH "Path to the Boost installation (or to the Boost source)")
option (VLI_DEBUG_COMPILE "Create a debug compile of VLI" OFF)
option (VLI_USE_GPU "Compile the GPU enhanced VLI components" ON)
option (VLI_USE_GMP "Compile with GMP components" ON)
option (VLI_TESTS "Build and execute the VLI regression tests" ON)
option (VLI_USE_OPENMP "Build and execute with openmp" ON)

set (VLI_CPU_FREQ "2.66" CACHE STRING "Frequency of the CPU (for benchmark timers)")
set (VLI_COMPILE_FLAGS_RELEASE " -Wall -O3 -funroll-loops -mtune=native -m64" CACHE STRING "Compiler flags for a regular compile.")
set (VLI_COMPILE_FLAGS_DEBUG " -Wall -g -O0 -m64" CACHE STRING "Compiler flags for a debug compile.")
set (VLI_GPU_TYPES "((unsigned int,3)) ((unsigned int,4)) ((unsigned long int,3)) ((unsigned long int,4))" CACHE STRING "What sizes and BaseInt types should be used for the vli_gpu class?")

#add_definitions(-DPOLYNOMIAL_MAX_ORDER=8 -DVLI_MAX_ORDER=4 -DVECTOR_MAX_ORDER=512)
add_definitions(-DUSE_VLI_INTEGERS_GPU -DPOLYNOMIAL_MAX_ORDER=21)

########################################################################
#
# Find dependencies
#
########################################################################

if(VLU_USE_OPENMP)
    add_definitions(-fopenmp)
    set(ENV{OMP_NUM_THREADS} "2")
endif(VLU_USE_OPENMP)


if (VLI_USE_GPU)
   find_package (CUDA) 
    include_directories(${CUDA_INCLUDE_DIRS})
endif (VLI_USE_GPU)

if (VLI_USE_GMP)

find_path(GMP_INCLUDE_DIR gmp.h $ENV{HOME}/.local/include)
find_library(GMP_LIBRARY_CXX NAMES gmpxx PATH $ENV{HOME}/.local/include)
find_library(GMP_LIBRARY_C NAMES gmp PATH $ENV{HOME}/.local/include)

if (GMP_INCLUDE_DIR AND GMP_LIBRARY_C AND GMP_LIBRARY_CXX)
   set(GMP_FOUND true)
endif (GMP_INCLUDE_DIR AND GMP_LIBRARY_C AND GMP_LIBRARY_CXX)

if (GMP_FOUND)
   if (NOT GMP_FIND_QUIETLY)
      message(STATUS "Found GMP: ${GMP_LIBRARY_CXX} ${GMP_LIBRARY_C}")
   endif (NOT GMP_FIND_QUIETLY)
else (GMP_FOUND)
   if (GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find GMP")
   endif (GMP_FIND_REQUIRED)
endif (GMP_FOUND)
endif (VLI_USE_GMP)

set(BOOST_ROOT ${BOOST_DIR})
find_package (Boost COMPONENTS unit_test_framework)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost libraries not found. Please specify location using the BOOST_ROOT variable")
endif()

########################################################################
#
# Configure Header File Template
#
########################################################################
configure_file (
        "${PROJECT_SOURCE_DIR}/vli_utils/vli_config.h.in"
        "${PROJECT_BINARY_DIR}/vli_utils/vli_config.h"
        )
include_directories ("${PROJECT_BINARY_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}")

########################################################################
#
# Kernel Compile
#
########################################################################

if(VLI_USE_GPU)
    CUDA_ADD_LIBRARY(kernels_gpu detail/kernels_gpu.cu)
endif(VLI_USE_GPU)

########################################################################
#
# Main Compile
#
########################################################################

add_executable (vli main.cpp)

if(VLI_USE_GPU)
    target_link_libraries (vli kernels_gpu ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES} ${GMP_LIBRARIES})
endif(VLI_USE_GPU)



########################################################################
#
# Tests
#
########################################################################

if(VLI_TESTS)
    enable_testing ()
    add_subdirectory (regression/)
endif(VLI_TESTS)


########################################################################
#
# Set internal variables according to the variables in section "options"
#
########################################################################
mark_as_advanced(CMAKE_BUILD_TYPE)
if(VLI_DEBUG_COMPILE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "CMake Build Type (Release/Debug)" FORCE)
else(VLI_DEBUG_COMPILE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake Build Type (Release/Debug)" FORCE)
endif(VLI_DEBUG_COMPILE)
set(CMAKE_CXX_FLAGS_RELEASE ${VLI_COMPILE_FLAGS_RELEASE} CACHE STRING "Compiler flags for Release compiles." FORCE)
set(CMAKE_CXX_FLAGS_DEBUG ${VLI_COMPILE_FLAGS_DEBUG} CACHE STRING "Compiler flags for Debug compiles." FORCE)

