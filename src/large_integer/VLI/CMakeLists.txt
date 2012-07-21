cmake_minimum_required (VERSION 2.8)

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
option (VLI_COMPILE_GPU "Compile the GPU enhanced VLI components" ON)
option (VLI_TESTS "Build the VLI regression tests" OFF)
option (VLI_BENCHMARKS "Build the VLI benchmarks" OFF)
option (VLI_USE_OPENMP "Build and execute with openmp" ON)
option (VLI_MAIN  "Build and execute the main" ON)
option (VLI_ARCH_X86  "Build VLI on x86 architecture  " ON)
option (VLI_ARCH_POWER  "Build VLI on Power architecture " OFF)

if(NOT MACHINE_CONFIG)
    set (VLI_COMPILER_FLAGS_RELEASE " -Wall -O2 -funroll-loops -m64" CACHE STRING "Compiler flags for a regular compile.")
    set (VLI_COMPILER_FLAGS_DEBUG " -Wall -g -O0 -m64" CACHE STRING "Compiler flags for a debug compile.")
endif(NOT MACHINE_CONFIG)

set (VLI_CPU_FREQ "2.5e9" CACHE STRING "Frequency of the CPU (for benchmark timers)")
set (VLI_SPLIT_PARAM "0.8" CACHE STRING "Split param for cpu/gpu mix mode, 1.0 = full gpu")
set (VLI_COMPILE_TYPES "((unsigned long int,3,10,3))" CACHE STRING "What sizes and BaseInt types should be compiled? Affects only the vli_gpu class and the tests.")

########################################################################
#
# Find dependencies
#
########################################################################
if(VLI_USE_OPENMP)
     find_package(OpenMP)
     add_definitions(${OpenMP_CXX_FLAGS})
endif(VLI_USE_OPENMP)

if (VLI_COMPILE_GPU)
   find_package (CUDA)
# set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch=sm_20")
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch=sm_20 --ptxas-options=-v") 
   include_directories(${CUDA_INCLUDE_DIRS})
else (VLI_COMPILE_GPU)
   set (CUDA_LIBRARIES "")
endif (VLI_COMPILE_GPU)

# The tests need the GMP library
#if (VLI_TESTS)
    find_path(GMP_INCLUDE_DIR gmp.h $ENV{HOME}/.local/include)
    find_library(GMP_LIBRARY_CXX NAMES gmpxx PATH $ENV{HOME}/.local/include)
    find_library(GMP_LIBRARY_C NAMES gmp PATH $ENV{HOME}/.local/include)

    if (GMP_INCLUDE_DIR AND GMP_LIBRARY_C AND GMP_LIBRARY_CXX)
       set(GMP_LIBRARY_CXX ${GMP_LIBRARY_CXX} ${GMP_LIBRARY_C})
       set(GMP_FOUND true)
       include_directories(${GMP_INCLUDE_DIR})
    endif (GMP_INCLUDE_DIR AND GMP_LIBRARY_C AND GMP_LIBRARY_CXX)

    if (GMP_FOUND)
      message(STATUS "Found GMP: ${GMP_LIBRARY_CXX} ${GMP_LIBRARY_C}")
    else (GMP_FOUND)
      message(WARNING "Could not find GMP library which are needed for the regression tests -> Deactivating tests.")
      set (VLI_TESTS OFF)
    endif (GMP_FOUND)
#endif (VLI_TESTS)

set(BOOST_ROOT ${BOOST_DIR})
find_package (Boost COMPONENTS unit_test_framework random)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost libraries not found. Please specify location using the BOOST_DIR variable")
endif()

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
set(CMAKE_CXX_FLAGS_RELEASE ${VLI_COMPILER_FLAGS_RELEASE} CACHE STRING "Compiler flags for Release compiles." FORCE)
set(CMAKE_CXX_FLAGS_DEBUG ${VLI_COMPILER_FLAGS_DEBUG} CACHE STRING "Compiler flags for Debug compiles." FORCE)


########################################################################
#
# Compile and prepare library
#
########################################################################
include_directories ("${PROJECT_BINARY_DIR}") 	
include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/src/")

add_subdirectory( src/vli/ )


########################################################################
#
# Main Compile
#
########################################################################
if(VLI_MAIN)
    add_subdirectory(main/)
endif()

########################################################################
#
# Tests
#
########################################################################

if(VLI_TESTS)
    include(CTest)
    enable_testing() 
    set(BUILDNAME "${PROJECT_VERSION}" CACHE STRING "Name of build on the dashboard")
      mark_as_advanced(BUILDNAME)
    add_subdirectory (regression)
    add_subdirectory (regression/vli_cpu)
endif(VLI_TESTS)

if(VLI_BENCHMARKS)
    include(CTest)
    enable_testing() 
    set(BUILDNAME "${PROJECT_VERSION}" CACHE STRING "Name of build on the dashboard")
      mark_as_advanced(BUILDNAME)
    add_subdirectory (benchmarks)
endif(VLI_BENCHMARKS)


########################################################################
#
# Install
#
########################################################################

install(DIRECTORY src/vli DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN ".svn" EXCLUDE)
install(FILES ${PROJECT_BINARY_DIR}/src/vli/libvli_cpu.a DESTINATION lib)
install(FILES ${PROJECT_BINARY_DIR}/vli/vli_config.h DESTINATION include/vli/)
if(VLI_COMPILE_GPU)
install(FILES ${PROJECT_BINARY_DIR}/src/vli/libvli_gpu.a DESTINATION lib)
endif(VLI_COMPILE_GPU)
