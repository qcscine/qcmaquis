cmake_minimum_required (VERSION 2.6)


########################################################################
#
# DEBUG mode
#
########################################################################
set(CMAKE_VERBOSE_MAKEFILE false)

project (vli)

########################################################################
#
# Version information
#
########################################################################

set (VLI_VERSION_MAJOR 0)
set (VLI_VERSION_MINOR 1)



########################################################################
#
# Options
#
########################################################################

set (BOOST_ROOT ${BOOST_ROOT} CACHE PATH "Path to the Boost installation (or to the Boost source)")
option (VLI_USE_GPU "Compile the GPU enhanced VLI components" ON)
option (VLI_TESTS "Build and execute the VLI regression tests" ON)
set (VLI_CPU_FREQ "2.66" CACHE STRING "Frequency of the CPU (for benchmark timers)")



########################################################################
#
# Find dependencies
#
########################################################################

if (VLI_USE_GPU)
    find_package (CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
endif (VLI_USE_GPU)

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
    target_link_libraries (vli kernels_gpu ${CUDA_CUBLAS_LIBRARIES}  ${CUDA_LIBRARIES})
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


