
cmake_minimum_required(VERSION 2.8)

######################################################################
# Machine Config
######################################################################
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
project (MAQUIS_TYPES) 

set(AMBIENT_SOURCES_DIR "../ambient" CACHE STRINGS "Location of the ambient source tree")

option(ENABLE_MTM "Enable multi threads matrix class" OFF)
option(USE_GPU "Use GPU in libmtm" OFF)
option(ENABLE_PARALLEL "Enable parallel matrix classes" OFF)
option(ENABLE_REGRESSION_VALIDATION "Build regression validation tests" OFF)
option(ENABLE_REGRESSION_PERFORMANCE "Build regression performance tests" OFF)
option(ENABLE_ALGEBRA "Enable algebra test " OFF)
option(ENABLE_PLOT "Plot performance test " OFF)
option(BUILD_AMBIENT "Build Ambient ? " OFF)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/maquis" CACHE PATH "MAQUIS install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


######################################################################
# Version information
######################################################################
set(PROJECT_YEAR 2011)
set(PROJECT_VERSION_MAJOR 2)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_BUILD "")

if(NOT PROJECT_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
  FIND_PACKAGE(Subversion) 
  IF(Subversion_FOUND) 
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} PROJECT) 
    set(PROJECT_VERSION_BUILD "r${PROJECT_WC_REVISION}")
  ENDIF(Subversion_FOUND) 
endif(NOT PROJECT_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

if(PROJECT_VERSION_BUILD)
  set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}-${PROJECT_VERSION_BUILD}")
else(PROJECT_VERSION_BUILD)
  set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
endif(PROJECT_VERSION_BUILD)
set(PROJECT_VERSION_STRING "Maquis - Matrices version ${PROJECT_VERSION}")
MESSAGE(STATUS "Maquis - Matrices version: ${PROJECT_VERSION}") 


######################################################################
# Libraries
######################################################################

#gnuplot
if(ENABLE_PLOT)
    find_package(Gnuplot)
endif(ENABLE_PLOT)

# Boost
set(Boost_requirements bindings)
if(BUILD_MTM)
  list(APPEND Boost_requirements thread)
endif(BUILD_MTM)
if(ENABLE_REGRESSION_VALIDATION OR ENABLE_REGRESSION_PERFORMANCE)
  list(APPEND Boost_requirements unit_test_framework)
endif(ENABLE_REGRESSION_VALIDATION OR ENABLE_REGRESSION_PERFORMANCE)

if(ALPS_ROOT_DIR)
  set(ALPS_DIR ${ALPS_ROOT_DIR}/share/alps)
  set(PREVENT_ALPS_COMPILERS ON)
  find_package(ALPS NO_SYSTEM_ENVIRONMENT_PATH)
endif(ALPS_ROOT_DIR)

include(${PROJECT_SOURCE_DIR}/config/find_boost.cmake)
#find_maquis_boost(COMPONENTS date_time filesystem program_options python regex system serialization)
find_maquis_boost(COMPONENTS ${Boost_requirements})

message(STATUS "Boost include: ${MAQUIS_Boost_INCLUDE_DIRS}")
message(STATUS "Boost libs: ${MAQUIS_Boost_LIBRARIES}")
message(STATUS "Boost unit_test: ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")

include_directories(${MAQUIS_Boost_INCLUDE_DIRS})
list(APPEND MAQUIS_TYPES_INCLUDE_DIRS ${MAQUIS_Boost_INCLUDE_DIRS})

link_directories(${MAQUIS_Boost_LIBRARY_DIRS})
list(APPEND MAQUIS_TYPES_LIBRARY_DIRS ${MAQUIS_Boost_LIBRARY_DIRS})

list(APPEND MAQUIS_TYPES_LIBRARIES ${MAQUIS_Boost_LIBRARIES})

# Blas/Lapack
if(ENABLE_ALGEBRA)
    if(!APPLE) # are you on a mac ?
        enable_language(Fortran) #needeed for blas/lapack
    endif(!APPLE)

    find_package(BLAS)
    find_package(LAPACK)
endif(ENABLE_ALGEBRA)

if(BUILD_MTM)
  # pthreads
  find_package(Threads REQUIRED)
  list(APPEND MAQUIS_TYPES_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
  
  if(USE_GPU)
    # Cuda
    find_package(CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
    list(APPEND MAQUIS_TYPES_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS})
    list(APPEND MAQUIS_TYPES_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
    add_definitions(-DUSE_GPU)
  endif(USE_GPU)
endif(BUILD_MTM)


if(ENABLE_PARALLEL)
  # first we need BLAS/LAPAACK
  find_package(MPI) 
  # second we need Ambient
 
  if(BUILD_AMBIENT)
    add_subdirectory(${AMBIENT_SOURCES_DIR} ambient-build)
    include(${CMAKE_BINARY_DIR}/ambient-build/MaquisAmbientBuildTree.cmake)
  else(BUILD_AMBIENT)
    if(NOT MaquisAmbient_FOUND)
      find_package(MaquisAmbient REQUIRED)
    endif(NOT MaquisAmbient_FOUND)
  endif(BUILD_AMBIENT)
  macro(compile_with_ambient)
    include_directories(${MAQUIS_AMBIENT_INCLUDE_DIRS})
    link_directories(${MAQUIS_AMBIENT_LIBRARY_DIRS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPI_PARALLEL ${MAQUIS_AMBIENT_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MAQUIS_AMBIENT_EXE_FLAGS}")
  endmacro(compile_with_ambient)
  else(USE_AMBIENT)
  macro(compile_with_ambient)
    message(ERROR Ambient not enabled!)
  endmacro(compile_with_ambient)

  #find ambient
  #include_directories(${MAQUIS_AMBIENT_INCLUDE_DIRS})
  #list(APPEND MAQUIS_TYPES_INCLUDE_DIRS ${MAQUIS_AMBIENT_INCLUDE_DIRS})
  #list(APPEND MAQUIS_TYPES_LIBRARY_DIRS ${MAQUIS_AMBIENT_LIBRARY_DIRS})
  #list(APPEND MAQUIS_TYPES_LIBRARIES ${MAQUIS_AMBIENT_LIBRARIES})

endif(ENABLE_PARALLEL)

######################################################################
# CPU Freq.
######################################################################
include(${PROJECT_SOURCE_DIR}/config/cpu_freq.cmake)


######################################################################
# How to build
######################################################################

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
  # Uncomment if you want the '-g' option in release mode (by default)
  # set(CMAKE_CXX_FLAGS_RELEASE "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "CMAKE_CXX_FLAGS_RELEASE: Flags used by the compiler during release builds" FORCE)
  # set(CMAKE_C_FLAGS_RELEASE "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "CMAKE_C_FLAGS_RELEASE: Flags used by the compiler during release builds" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})



######################################################################
# Targets
######################################################################

# Include/Link Directories
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)
list(APPEND MAQUIS_TYPES_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include)
list(APPEND MAQUIS_TYPES_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib)


# Targets
add_subdirectory(src)


# ALPS I need ALPS for Boost numeric
set(ALPS_DIR $ENV{ALPS_ROOT} CACHE PATH "Path to the Alps include, need for boost numeric")

if(ENABLE_REGRESSION_VALIDATION OR ENABLE_REGRESSION_PERFORMANCE)
  include(CTest)
  enable_testing() 
endif(ENABLE_REGRESSION_VALIDATION OR ENABLE_REGRESSION_PERFORMANCE)

if(ENABLE_REGRESSION_VALIDATION)
  add_subdirectory(regression/dense_matrix)
  if(ENABLE_PARALLEL)
    add_subdirectory(regression/validation)
  endif(ENABLE_PARALLEL)
  if(ENABLE_MTM)
    add_subdirectory(regression/mt_matrix)
  endif(ENABLE_MTM)
endif(ENABLE_REGRESSION_VALIDATION)

if(ENABLE_REGRESSION_PERFORMANCE)
  if(ENABLE_PARALLEL)
    add_subdirectory(regression/performance)
  endif(ENABLE_PARALLEL)
endif(ENABLE_REGRESSION_PERFORMANCE)

# Configure
configure_file(MaquisTypesConfig.cmake.in ${CMAKE_BINARY_DIR}/MaquisTypesConfig.cmake)


# Install / export
install(FILES ${CMAKE_BINARY_DIR}/MaquisTypesConfig.cmake
        DESTINATION share COMPONENT dev)
if(ENABLE_MTM)
  install(EXPORT MaquisTypesTargets  DESTINATION share COMPONENT dev)
endif(ENABLE_MTM)
