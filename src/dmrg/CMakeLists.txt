
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
project (MAQUIS_DMRG) 

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/maquis" CACHE PATH "MAQUIS install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(MATRICES_DIR "../types/src" CACHE STRINGS "Location of the matrices sources")
set(AMBIENT_SOURCES_DIR "../ambient" CACHE STRINGS "Location of the ambient source tree")

set(BUILD_SYMMETRIES "U1;TwoU1;NONE" CACHE STRINGS "List of symmetry to include in the build objects. U1, TwoU1, NONE are currently available.")

option(ENABLE_ALPS_MODELS "Enable binding with ALPS lattices and models" ON)
option(ENABLE_OMP "Enable OpenMP" OFF)
option(USE_AMBIENT "Enable ambient" OFF)
option(USE_MTM "Use multi threads matrix class" OFF)
option(USE_GPU "Use GPU in libmtm" OFF)
option(BUILD_AMBIENT "Build Ambient inside this build tree" OFF)

option(ENABLE_COLLECTOR "Enable profiling through DataCollector object" OFF)

option(BUILD_REGRESSION "Build regression tests" OFF)

mark_as_advanced(ENABLE_COLLECTOR)

set(BLAS_LAPACK_SELECTOR ${DEFAULT_BLAS_LAPACK} CACHE STRING "Choose the Blas/Lapack version to use, possible values 'alps', 'veclib', 'manual'")
set(BLAS_LAPACK_MANUAL_LIBS ${BLAS_LAPACK_MANUAL_LIBS} CACHE STRINGS "")
set(BLAS_LAPACK_MANUAL_LIBS_DIR ${BLAS_LAPACK_MANUAL_LIBS_DIR} CACHE STRINGS "")
set(BLAS_LAPACK_MANUAL_DEFINITIONS ${BLAS_LAPACK_MANUAL_DEFINITIONS} CACHE STRINGS "")
set(BLAS_LAPACK_MANUAL_INCLUDES ${BLAS_LAPACK_MANUAL_INCLUDES} CACHE STRINGS "")


######################################################################
# How to build
######################################################################

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "CMAKE_CXX_FLAGS_RELEASE: Flags used by the compiler during release builds" FORCE)
  set(CMAKE_C_FLAGS_RELEASE " ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "CMAKE_C_FLAGS_RELEASE: Flags used by the compiler during release builds" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


######################################################################
# Libraries
######################################################################
# OpenMP
if(NOT OpenMP_FOUND)
  find_package(OpenMP)
endif(NOT OpenMP_FOUND)

# ALPS
if(NOT ALPS_DIR AND ALPS_ROOT_DIR)
  set(ALPS_DIR ${ALPS_ROOT_DIR}/share/alps)
endif(NOT ALPS_DIR AND ALPS_ROOT_DIR)

set(PREVENT_ALPS_COMPILERS ON)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)

set(ALPS_LIBRARY_DIRS ${ALPS_LIBRARY_DIRS} CACHE STRINGS "ALPS library dirs")
set(ALPS_INCLUDE_DIRS ${ALPS_INCLUDE_DIRS} CACHE STRINGS "ALPS include dirs")
if(ALPS_HDF5_INCLUDE_DIR) # support for old ALPS revisions
  set(HDF5_INCLUDE_DIR ${ALPS_HDF5_INCLUDE_DIR} CACHE STRINGS "")
endif(ALPS_HDF5_INCLUDE_DIR)
if(ALPS_HDF5_LIBRARIES) # support for old ALPS revisions
  set(HDF5_LIBRARIES ${ALPS_HDF5_LIBRARIES} CACHE STRINGS "")
endif(ALPS_HDF5_LIBRARIES)

if (LAPACK_FOUND) # support for old ALPS revisions
  set(ALPS_LAPACK_DEFINITIONS ${LAPACK_DEFINITIONS} CACHE STRINGS "")
  set(ALPS_LAPACK_LINKER_FLAGS "${LAPACK_LINKER_FLAGS}" CACHE STRINGS "")
  set(ALPS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRINGS "")
  set(ALPS_LAPACK_LIBRARY ${LAPACK_LIBRARY} CACHE STRINGS "")
  set(ALPS_BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE STRINGS "")
  set(ALPS_BLAS_LIBRARY ${BLAS_LIBRARY} CACHE STRINGS "")
endif(LAPACK_FOUND)

mark_as_advanced(ALPS_LIBRARY_DIRS ALPS_INCLUDE_DIRS HDF5_INCLUDE_DIR HDF5_LIBRARIES)


# Boost
include(${PROJECT_SOURCE_DIR}/config/boost.cmake)
include_directories(${MAQUIS_Boost_INCLUDE_DIRS})
link_directories(${MAQUIS_Boost_LIBRARY_DIRS})

# Ambient
if(USE_AMBIENT)
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
endif(USE_AMBIENT)

######################################################################
# Version information
######################################################################
set(DMRG_YEAR 2011)
set(DMRG_VERSION_MAJOR 2)
set(DMRG_VERSION_MINOR 0)
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
include(${PROJECT_SOURCE_DIR}/config/cpu_freq.cmake)


######################################################################
# Parsing user options
######################################################################

set(DMRG_LIBRARIES alps ${MAQUIS_Boost_LIBRARIES} ${HDF5_LIBRARIES})


# BLAS/LAPACK selector
if(${BLAS_LAPACK_SELECTOR} MATCHES "alps")
  if(CMAKE_EXE_LINKER_FLAGS MATCHES "^ $")
    set(CMAKE_EXE_LINKER_FLAGS ${ALPS_LAPACK_LINKER_FLAGS})
  else(CMAKE_EXE_LINKER_FLAGS MATCHES "^ $")
    list(APPEND CMAKE_EXE_LINKER_FLAGS ${ALPS_LAPACK_LINKER_FLAGS})
  endif(CMAKE_EXE_LINKER_FLAGS MATCHES "^ $")
  set(BLAS_LAPACK_LIBRARIES ${ALPS_LAPACK_LIBRARIES} ${ALPS_LAPACK_LIBRARY}
      ${ALPS_BLAS_LIBRARIES} ${ALPS_BLAS_LIBRARY})
  list(APPEND DMRG_LIBRARIES ${BLAS_LAPACK_LIBRARIES})
  add_definitions(${ALPS_LAPACK_DEFINITIONS})
endif(${BLAS_LAPACK_SELECTOR} MATCHES "alps")
if(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework vecLib")
    set(BLAS_LAPACK_LIBRARIES /usr/lib/libpthread.dylib)
    list(APPEND DMRG_LIBRARIES ${BLAS_LAPACK_LIBRARIES})
  else(APPLE)
    meassage("Error: VecLib available only on Mac.")
  endif(APPLE)
endif(${BLAS_LAPACK_SELECTOR} MATCHES "veclib")
if(${BLAS_LAPACK_SELECTOR} MATCHES "manual")
  set(BLAS_LAPACK_LIBRARIES ${BLAS_LAPACK_MANUAL_LIBS})
  list(APPEND DMRG_LIBRARIES ${BLAS_LAPACK_LIBRARIES})
  add_definitions(${BLAS_LAPACK_MANUAL_DEFINITIONS})
  include_directories(${BLAS_LAPACK_MANUAL_INCLUDES})
  link_directories(${BLAS_LAPACK_MANUAL_LIBS_DIR})
endif(${BLAS_LAPACK_SELECTOR} MATCHES "manual")

# LIST(REVERSE CMAKE_EXE_LINKER_FLAGS)
# LIST(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
# LIST(REVERSE CMAKE_EXE_LINKER_FLAGS)
# message("CMAKE_EXE_LINKER_FLAGS:"${CMAKE_EXE_LINKER_FLAGS})


if(USE_MTM)
  add_definitions(-DUSE_MTM)
  if(USE_GPU)
    find_package(CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
    list(APPEND DMRG_LIBRARIES ${CUDA_CUBLAS_LIBRARIES} ${CUDA_LIBRARIES})
    add_definitions(-DUSE_GPU)
  endif(USE_GPU)
else(USE_MTM)
  if(USE_GPU)
    message(STATUS "USE_GPU only works if USE_MTM is enabled")
  endif(USE_GPU)
endif(USE_MTM)



# OpenMP
if(ENABLE_OMP)
  if(OPENMP_FOUND)
    macro(enable_omp_if_found)
      #add_definitions(${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endmacro(enable_omp_if_found)
  else(OPENMP_FOUND)
    message("OpenMP not found")
  endif(OPENMP_FOUND)
else(ENABLE_OMP)
  macro(enable_omp_if_found)
  endmacro(enable_omp_if_found)
endif(ENABLE_OMP)


######################################################################
# Symmetries management
######################################################################

macro(get_symm_suffix RET SYMM)
  set(${RET} ${SYMM})
  string(REGEX REPLACE "^U1$" "u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1$" "2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1$" "2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NONE$" "none" ${RET} ${${RET}})
endmacro(get_symm_suffix)

macro(get_symm_group_name RET SYMM)
  set(${RET} ${SYMM})
  string(REGEX REPLACE "^U1$" "U1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1$" "TwoU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1$" "TwoU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NONE$" "TrivialGroup" ${RET} ${${RET}})
endmacro(get_symm_group_name)

macro(get_symm_files TYPE RET FILEBASE)
  if(NOT ${TYPE} STREQUAL "APPEND")
    set(${RET} "")
  endif(NOT ${TYPE} STREQUAL "APPEND")
  foreach(SYMM ${BUILD_SYMMETRIES})
    get_symm_suffix(SYMM_SUFFIX ${SYMM})
    string(REPLACE "{SYMM}" ${SYMM_SUFFIX} SYMM_FILE ${FILEBASE})
    list(APPEND ${RET} ${SYMM_FILE})
  endforeach(SYMM)
endmacro(get_symm_files)

macro(configure_symm_file INPUT OUTBASE VARNAME)
  foreach(SYMM ${BUILD_SYMMETRIES})
    get_symm_suffix(SYMM_SUFFIX ${SYMM})
    get_symm_group_name(${VARNAME} ${SYMM})
    string(REPLACE "{SYMM}" ${SYMM_SUFFIX} SYMM_FILE ${OUTBASE})
    configure_file(${INPUT} ${SYMM_FILE})
  endforeach(SYMM)
endmacro(configure_symm_file)


# Adding symmetries definitions
foreach(SYMM ${BUILD_SYMMETRIES})
  get_symm_group_name(SYMM_NAME ${SYMM})
  add_definitions(-DHAVE_${SYMM_NAME})
endforeach(SYMM)




######################################################################
# Include / link directories
######################################################################

include_directories(${PROJECT_SOURCE_DIR}/framework ${PROJECT_BINARY_DIR}/framework
                    ${MATRICES_DIR}
                    ${HDF5_INCLUDE_DIR} ${ALPS_INCLUDE_DIRS})
link_directories(${ALPS_LIBRARY_DIRS})
set(DMRG_FRAMEWORK_DIR ${PROJECT_SOURCE_DIR}/framework)


######################################################################
# Configure files
######################################################################

configure_file(framework/dmrg_version.h.in ${CMAKE_BINARY_DIR}/framework/dmrg_version.h)


######################################################################
# Targets
######################################################################

# DataCollector
if(ENABLE_COLLECTOR)
  add_definitions(-DENABLE_DATACOLLECTORS)
  add_library(maquis_collector framework/utils/data_collector.cpp)
  list(APPEND DMRG_LIBRARIES maquis_collector)
endif(ENABLE_COLLECTOR)

# *** libraries / utilities
if(USE_MTM)
add_subdirectory(${MATRICES_DIR}/types/mt_matrix types/mt_matrix)
endif(USE_MTM)
add_subdirectory(framework/dmrg)

# *** applications
add_subdirectory(applications)

# *** regression tests
if(BUILD_REGRESSION)
  enable_testing()
	if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
		message(WARNING "Tests should be checked in Debug mode!")
	endif()
  add_subdirectory(regression/block_matrix)
  add_subdirectory(regression/mp_tensors)
  add_subdirectory(regression/dmrg)
  add_subdirectory(regression/te)
  if(USE_AMBIENT)
     add_subdirectory(regression/pdmrg)
  endif(USE_AMBIENT)
else(BUILD_REGRESSION)
  add_subdirectory(regression/dmrg EXCLUDE_FROM_ALL)
  if(USE_AMBIENT)
     add_subdirectory(regression/pdmrg EXCLUDE_FROM_ALL)
  endif(USE_AMBIENT)
endif(BUILD_REGRESSION)



######################################################################
# Install
######################################################################


install(FILES ${CMAKE_BINARY_DIR}/framework/dmrg_version.h
        DESTINATION include COMPONENT headers)

install(DIRECTORY framework/dmrg framework/ietl framework/utils
        DESTINATION include
        COMPONENT headers
        PATTERN ".svn" EXCLUDE
        PATTERN ".DS_Store" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.C" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.cu" EXCLUDE
        PATTERN "*.in" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN "Makefile*" EXCLUDE
        )


