cmake_minimum_required (VERSION 2.6)

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

set (BOOST_ROOT "/opt/boost" CACHE PATH "Path to the Boost installation (or to the Boost source)")
option (VLI_USE_GPU "Compile the GPU enhanced VLI components" ON)
option (VLI_TESTS "Build and execute the VLI regression tests" ON)
set (VLI_CPU_FREQ "2.66" CACHE STRING "Frequency of the CPU (for benchmark timers)")



########################################################################
#
# Find dependencies
#
########################################################################

# include our cmake modules if they are not present in cmake installation
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-files)

if (VLI_USE_GPU)
    find_package (CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
endif (VLI_USE_GPU)



########################################################################
#
# Configure Header File Template
#
########################################################################
configure_file (
        "${PROJECT_SOURCE_DIR}/vli_config.h.in"
        "${PROJECT_BINARY_DIR}/vli_config.h"
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


