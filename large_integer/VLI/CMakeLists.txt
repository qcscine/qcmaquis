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
option (VLI_TESTS "Build the VLI regression tests" ON)
option (VLI_USE_OPENMP "Build and execute with openmp" ON)

set (VLI_CPU_FREQ "2.5e9" CACHE STRING "Frequency of the CPU (for benchmark timers)")
set (VLI_SPLIT_PARAM "0.1" CACHE STRING "Split param for cpu/gpu mix mode ")
set (VLI_COMPILER_FLAGS_RELEASE " -Wall -O2 -funroll-loops -m64" CACHE STRING "Compiler flags for a regular compile.")
set (VLI_COMPILER_FLAGS_DEBUG " -Wall -g -O0 -m64" CACHE STRING "Compiler flags for a debug compile.")
set (VLI_COMPILE_TYPES "((unsigned int,3)) ((unsigned int,4)) ((unsigned long int,3)) ((unsigned long int,4))" CACHE STRING "What sizes and BaseInt types should be compiled? Affects only the vli_gpu class and the tests.")

add_definitions(-DUSE_VLI_INTEGERS_GPU -DPOLYNOMIAL_MAX_ORDER=21)

########################################################################
#
# Find dependencies
#
########################################################################
if(VLI_USE_OPENMP)
    add_definitions(-fopenmp)
# TODO is there a more general way to find the right library?
    set(OPENMP_LIBRARIES "gomp")
endif(VLI_USE_OPENMP)

if (VLI_USE_GPU)
   find_package (CUDA)
   include_directories(${CUDA_INCLUDE_DIRS})
else (VLI_USE_GPU)
   set (CUDA_LIBRARIES "")
endif (VLI_USE_GPU)

# The tests need the GMP library
if (VLI_TESTS)
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
endif (VLI_TESTS)

set(BOOST_ROOT ${BOOST_DIR})
find_package (Boost COMPONENTS unit_test_framework)
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

add_executable (a.out main.cpp)

if(VLI_USE_GPU)
    target_link_libraries (a.out vli_gpu ${CUDA_LIBRARIES} ${GMP_LIBRARY_CXX} ${OPENMP_LIBRARIES})
else(VLI_USE_GPU)
    target_link_libraries (a.out ${GMP_LIBRARY_CXX} ${OPENMP_LIBRARIES})
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
# Install
#
########################################################################

install(DIRECTORY src/vli DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN ".svn" EXCLUDE)
install(FILES ${PROJECT_BINARY_DIR}/vli/vli_config.h DESTINATION include/vli/)
install(FILES ${PROJECT_BINARY_DIR}/src/vli/libvli_gpu.a DESTINATION lib)

