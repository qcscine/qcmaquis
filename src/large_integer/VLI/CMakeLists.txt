cmake_minimum_required (VERSION 2.6)


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
option (VLI_TESTS "Build and execute the VLI regression tests" ON)
set (VLI_CPU_FREQ "2.66" CACHE STRING "Frequency of the CPU (for benchmark timers)")
set (VLI_COMPILE_FLAGS_RELEASE "-frecord-gcc-switches -Wall -O3 -funroll-loops -march=native" CACHE STRING "Compiler flags for a regular compile.")
set (VLI_COMPILE_FLAGS_DEBUG "-frecord-gcc-switches -Wall -g -O0" CACHE STRING "Compiler flags for a debug compile.")


########################################################################
#
# Find dependencies
#
########################################################################

if (VLI_USE_GPU)
    find_package (CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
endif (VLI_USE_GPU)

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
    target_link_libraries (vli kernels_gpu ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES} ${GMP_LIBRARIES} -lgmp)
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

