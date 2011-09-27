cmake_minimum_required (VERSION 2.6)

########################################################################
#
# Project and version information
#
########################################################################

project (se)
set (SE_VERSION_MAJOR 0)
set (SE_VERSION_MINOR 1)


########################################################################
#
# Options
#
########################################################################

set (BOOST_DIR $ENV{BOOST_DIR} CACHE PATH "Path to the Boost installation (or to the Boost source)")
set (VLI_DIR CACHE PATH "VLI install path")
 
option (SE_DEBUG_COMPILE "Create a debug compile of SE" OFF)
option (SE_USE_GMP " Compile the GMP version " ON)
option (SE_USE_VLI " Compile the VLI version " ON)
option (SE_USE_MIX " Compile the mix SE components" ON)
option (SE_USE_GPU " Compile the GPU enhanced SE components" ON)

set (SE_CPU_FREQ "2.5e9" CACHE STRING "Frequency of the CPU (for benchmark timers)")
set (SE_SPLIT_PARAM "0.1" CACHE STRING "Split param for cpu/gpu mix mode ")
set (SE_COMPILER_FLAGS_RELEASE " -Wall -O2 -funroll-loops -m64" CACHE STRING "Compiler flags for a regular compile.")
set (SE_COMPILER_FLAGS_DEBUG " -Wall -g -O0 -m64" CACHE STRING "Compiler flags for a debug compile.")

add_definitions(-DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG)

########################################################################
#
# Find dependencies
#
########################################################################

if (SE_USE_GPU)
   find_package (CUDA) 
   include_directories(${CUDA_INCLUDE_DIRS})
endif (SE_USE_GPU)

# The tests need the GMP library
if (SE_USE_GMP)
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
      set (SE_TESTS OFF)
    endif (GMP_FOUND)
endif (SE_USE_GMP)

#BOOST
set(BOOST_ROOT ${BOOST_DIR})
find_package (Boost COMPONENTS unit_test_framework)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost libraries not found. Please specify location using the BOOST_DIR variable")
endif()

#VLI
if(SE_USE_VLI)
    set(VLI_INCLUDE_DIR ${VLI_DIR}/include)  
    find_library(VLI_LIBRARY vli_gpu ${VLI_DIR}/lib)  

    if (VLI_INCLUDE_DIR AND VLI_LIBRARY)
       set(VLI_LIBRARY ${VLI_LIBRARY})
       set(VLI_FOUND true)
       include_directories(${VLI_INCLUDE_DIR})
    endif(VLI_INCLUDE_DIR AND VLI_LIBRARY)

    if(VLI_FOUND)
        message(STATUS "Found VLI library: ${VLI_LIBRARY}")
    else()
      message(WARNING "Could not find VLI library which is needed for the mix CPU/GPU version.")
    endif(VLI_FOUND)
endif(SE_USE_VLI)



########################################################################
#
# Compile and prepare library
#
########################################################################

include_directories ("${PROJECT_BINARY_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}")

########################################################################
#
# Main Compile
#
########################################################################

if(SE_USE_GMP)
   add_subdirectory(se_gmp)
endif(SE_USE_GMP)

if(SE_USE_VLI)
   add_subdirectory(se_vli)
endif(SE_USE_VLI)

if(SE_USE_MIX)
   add_subdirectory(se_vli_mix)
endif(SE_USE_MIX)

########################################################################
#
# Install
#
########################################################################

install(DIRECTORY src/vli DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN ".svn" EXCLUDE)
install(FILES ${PROJECT_BINARY_DIR}/vli/vli_config.h DESTINATION include/vli/)
install(FILES ${PROJECT_BINARY_DIR}/src/vli/libvli_gpu.a DESTINATION lib)

########################################################################
#
# Set internal variables according to the variables in section "options"
#
########################################################################

mark_as_advanced(CMAKE_BUILD_TYPE)
if(SE_DEBUG_COMPILE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "CMake Build Type (Release/Debug)" FORCE)
else(SE_DEBUG_COMPILE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake Build Type (Release/Debug)" FORCE)
endif(SE_DEBUG_COMPILE)
set(CMAKE_CXX_FLAGS_RELEASE ${SE_COMPILER_FLAGS_RELEASE} CACHE STRING "Compiler flags for Release compiles." FORCE)
set(CMAKE_CXX_FLAGS_DEBUG ${SE_COMPILER_FLAGS_DEBUG} CACHE STRING "Compiler flags for Debug compiles." FORCE)




