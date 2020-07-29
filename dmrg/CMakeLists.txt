cmake_minimum_required(VERSION 3.1.3)

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

if(NOT DEFAULT_BLAS_LAPACK)
    set(DEFAULT_BLAS_LAPACK auto)
endif(NOT DEFAULT_BLAS_LAPACK)

######################################################################
# Project START
######################################################################
project (MAQUIS_DMRG)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# set proper rpath for installed binaries
# @stknecht: make sure this works also for OS X with MACOS_RPATH
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt/maquis" CACHE PATH "MAQUIS install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set (ALPS_SOURCE_DIR "alps")

set(BUILD_SYMMETRIES "TwoU1;TwoU1PG;SU2U1;SU2U1PG" CACHE STRING "List of symmetry
to include in the build objects. U1, TwoU1, TwoU1PG, NU1, Z2, SU2U1, SU2U1PG and NONE are currently available.")
set(DMRG_NUMSYMM 6 CACHE INT "Maximum number of U1 symmetries.")
option(ENABLE_OMP "Enable OpenMP" ON)
# Disabled
# option(ENABLE_ALPS_MODELS "Enable binding with ALPS lattices and models" OFF)
option(LAPACK_64_BIT "Build for 64-bit LAPACK library" ON)

# Standalone utility binaries, mostly deprecated by the binary interface
option (BUILD_INIT "Build dmrg_init binary (deprecated)" OFF)
option (BUILD_CIDEAS "Build CIDEAS binary (deprecated)" OFF)
option (BUILD_MPS_ROTATE "Build mps_rotate binary for RASSI counterrotation (deprecated)" OFF)
option (BUILD_MPS_COMPRESS "Build mps_compress binary for MPS compression (not working currently)" OFF)
option (BUILD_MPS_OVERLAP "Build MPS Overlap calculator" OFF)
option (BUILD_ONETDM_DIAG "Build 1-TDM diagonal elements calculator (deprecated by dmrg_meas and the binary interface)" OFF)
option (BUILD_MPS2CI "Build CI coefficient reconstruction facilities" OFF)
option (BUILD_DEAS "Build s1 and I printer tool" OFF)
option (BUILD_MPS_TRANSFORM "Build MPS SU2->2U1 symmetry group transformation tool" OFF)

option(ENABLE_COLLECTOR "Enable profiling through DataCollector object" OFF)
set(BENCHMARKS_DIR "../../benchmarks" CACHE PATH "Location of benchmarks directory.")
option(TESTS      "Build and run tests" OFF)


mark_as_advanced(ENABLE_COLLECTOR)

set(BLAS_LAPACK_SELECTOR "auto" CACHE STRINGS "auto, manual, mkl_sequential, mkl_parallel, veclib, openblas")

######################################################################
# How to build
######################################################################
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  include("config/clang_fixes.cmake")
endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "-Wno-unknown-pragmas ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "CMAKE_CXX_FLAGS_RELEASE: Flags used by the compiler during release builds" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Wno-unknown-pragmas ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "CMAKE_CXX_FLAGS_RELWITHDEBINFO: Flags used by the compiler during Release with Debug Info builds." FORCE)
  set(CMAKE_C_FLAGS_RELEASE " ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "CMAKE_C_FLAGS_RELEASE: Flags used by the compiler during release builds" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


######################################################################
# Symmetries management
######################################################################

macro(get_symm_suffix RET SYMM)
  set(${RET} ${SYMM})
  string(REGEX REPLACE "^U1$" "u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1$" "2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1$" "2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NU1$"  "nu1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1PG$" "2u1pg" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1PG$" "2u1pg" ${RET} ${${RET}})
  string(REGEX REPLACE "^NONE$" "none" ${RET} ${${RET}})
  string(REGEX REPLACE "^Z2$" "Ztwo" ${RET} ${${RET}})
  string(REGEX REPLACE "^SU2U1$" "su2u1" ${RET} ${${RET}})
  string(REGEX REPLACE "^SU2U1PG$" "su2u1pg" ${RET} ${${RET}})
  string(REGEX REPLACE "^U1DG$" "u1dg" ${RET} ${${RET}})
endmacro(get_symm_suffix)

macro(get_symm_group_name RET SYMM)
  set(${RET} ${SYMM})
  string(REGEX REPLACE "^U1$" "U1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1$" "TwoU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1$" "TwoU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^NU1$"  "NU1" ${RET} ${${RET}})
  string(REGEX REPLACE "^TWOU1PG$" "TwoU1PG" ${RET} ${${RET}})
  string(REGEX REPLACE "^TwoU1PG$" "TwoU1PG" ${RET} ${${RET}})
  string(REGEX REPLACE "^NONE$" "TrivialGroup" ${RET} ${${RET}})
  string(REGEX REPLACE "^Z2$" "Ztwo" ${RET} ${${RET}})
  string(REGEX REPLACE "^SU2U1$" "SU2U1" ${RET} ${${RET}})
  string(REGEX REPLACE "^SU2U1PG$" "SU2U1PG" ${RET} ${${RET}})
  string(REGEX REPLACE "^U1DG$" "U1DG" ${RET} ${${RET}})
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
message(STATUS "Enabled symmetries: ${BUILD_SYMMETRIES}")
foreach(SYMM ${BUILD_SYMMETRIES})
  get_symm_group_name(SYMM_NAME ${SYMM})
  list(APPEND DMRG_DEFINITIONS -DHAVE_${SYMM_NAME})
  set(DMRG_HAS_${SYMM_NAME} TRUE)
endforeach(SYMM)

##### Dependencies #####

# GSL
# find GNU Scientific Library if using SU2 Groups
if((DMRG_HAS_SU2U1) OR (DMRG_HAS_SU2U1PG))
find_package(GSL REQUIRED)
  list(APPEND DMRG_LIBRARIES ${GSL_LIBRARIES})
  endif()

# BLAS/LAPACK selector
find_package(MAQUISLapack)
message(STATUS "Lapack include:  ${MAQUISLapack_INCLUDE_DIRS}")
message(STATUS "Lapack lib dirs: ${MAQUISLapack_LIB_DIRS}")
message(STATUS "Lapack libs:     ${MAQUISLapack_LIBRARIES}")
list(APPEND DMRG_LIBRARIES ${MAQUISLapack_LIBRARIES})
# ILP64 definition
if(LAPACK_64_BIT)
add_definitions(-DBIND_FORTRAN_INTEGER_8)
endif()

# HDF5
find_package(HDF5 REQUIRED)
list(APPEND DMRG_LIBRARIES ${HDF5_LIBRARIES})

# Boost
set(Boost_requirements program_options filesystem system serialization thread)
set (Boost_NO_BOOST_CMAKE ON)

find_package(Boost COMPONENTS ${Boost_requirements})
if(Boost_FOUND)
list(APPEND DMRG_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
list(APPEND DMRG_LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
list(APPEND DMRG_LIBRARIES ${Boost_LIBRARIES})
else(Boost_FOUND)
message(FATAL_ERROR "Boost libraries are required for QCMaquis")
endif(Boost_FOUND)

# ALPS
add_subdirectory(${ALPS_SOURCE_DIR})
list(APPEND DMRG_LIBRARIES alps)

######################################################################
# Version information
######################################################################
set(DMRG_YEAR 2018)
set(DMRG_VERSION_MAJOR 2)
set(DMRG_VERSION_MINOR 1)
set(DMRG_VERSION_BUILD "")

if(NOT DMRG_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn)
  FIND_PACKAGE(Subversion)
  IF(Subversion_FOUND)
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} DMRG)
    string(REPLACE ${DMRG_WC_ROOT} "" DMRG_BRANCH ${DMRG_WC_URL})
    set(DMRG_VERSION_BUILD "r${DMRG_WC_REVISION} (${DMRG_BRANCH})")
  ENDIF(Subversion_FOUND)
endif(NOT DMRG_VERSION_BUILD AND EXISTS ${PROJECT_SOURCE_DIR}/.svn)

if(EXISTS ${PROJECT_SOURCE_DIR}/../.git)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
            )
        if(NOT ${GIT_BRANCH} STREQUAL "")
            string(STRIP ${GIT_BRANCH} GIT_BRANCH)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE SVN_REV
            ERROR_QUIET
            )

        if(NOT ${SVN_REV} STREQUAL "")
            string(STRIP ${SVN_REV} SVN_REV)
        endif()

        set(DMRG_VERSION_BUILD "${GIT_REVISION} (${GIT_BRANCH}) ${SVN_REV}")
    endif()
endif(EXISTS ${PROJECT_SOURCE_DIR}/../.git)

if(DMRG_VERSION_BUILD)
  set(DMRG_VERSION "${DMRG_VERSION_MAJOR}.${DMRG_VERSION_MINOR}-${DMRG_VERSION_BUILD}")
else(DMRG_VERSION_BUILD)
  set(DMRG_VERSION "${DMRG_VERSION_MAJOR}.${DMRG_VERSION_MINOR}")
endif(DMRG_VERSION_BUILD)
set(DMRG_VERSION_STRING "DMRG version ${DMRG_VERSION}")
MESSAGE(STATUS "DMRG version: ${DMRG_VERSION}")

# OpenMP
if(ENABLE_OMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DMAQUIS_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else(OPENMP_FOUND)
    message("OpenMP not found")
  endif(OPENMP_FOUND)
endif(ENABLE_OMP)


######################################################################
# Include / link directories
######################################################################
set (ALPS_INCLUDE_DIRS ${ALPS_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/${ALPS_SOURCE_DIR}/src)
list(APPEND DMRG_INCLUDE_DIRS ${MATRICES_DIR} ${HDF5_INCLUDE_DIR} ${ALPS_INCLUDE_DIRS})
list(APPEND DMRG_LIBRARY_DIRS ${ALPS_LIBRARY_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/framework ${PROJECT_BINARY_DIR}/framework ${PROJECT_SOURCE_DIR}/lib/maquis_dmrg ${PROJECT_BINARY_DIR}/lib/maquis_dmrg ${DMRG_INCLUDE_DIRS})
link_directories(${DMRG_LIBRARY_DIRS})
add_definitions(${DMRG_DEFINITIONS})

set(DMRG_FRAMEWORK_DIR ${PROJECT_SOURCE_DIR}/framework)
if(CMAKE_EXE_LINKER_FLAGS MATCHES "^ $")
    set(CMAKE_EXE_LINKER_FLAGS ${MAQUISLapack_LINKER_FLAGS})
else(CMAKE_EXE_LINKER_FLAGS MATCHES "^ $")
    list(APPEND CMAKE_EXE_LINKER_FLAGS ${MAQUISLapack_LINKER_FLAGS})
endif(CMAKE_EXE_LINKER_FLAGS MATCHES "^ $")


######################################################################
# Configure files
######################################################################

configure_file(framework/dmrg/version.h.in ${CMAKE_BINARY_DIR}/framework/dmrg/version.h)


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
add_subdirectory(framework/dmrg)

# *** shared library
add_subdirectory(lib/maquis_dmrg)

# *** applications
add_subdirectory(applications)

# *** tests
if(TESTS)
  enable_testing()
  add_subdirectory(tests)
  add_test(NAME Test1 COMMAND test1)
  add_test(NAME Test2 COMMAND test2)
  add_test(NAME Test_Relativistic COMMAND test_rel)
  add_test(NAME Test_Integral_Map COMMAND test_integral_map)
  add_test(NAME Test_4rdm_save COMMAND test_4rdm_save)
endif(TESTS)




######################################################################
# DMRG exports
######################################################################

## For build tree
set(CONF_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/framework ${PROJECT_BINARY_DIR}/framework ${DMRG_INCLUDE_DIRS})
configure_file(MAQUIS_DMRGConfig.cmake.in ${PROJECT_BINARY_DIR}/MAQUIS_DMRGConfig.cmake @ONLY)

## For install
set(CONF_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include ${DMRG_INCLUDE_DIRS})
configure_file(MAQUIS_DMRGConfig.cmake.in ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/MAQUIS_DMRGConfig.cmake @ONLY)

install(FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/MAQUIS_DMRGConfig.cmake DESTINATION share/cmake  COMPONENT dev)
install(EXPORT DMRGTargets DESTINATION share/cmake COMPONENT dev)


######################################################################
# Install
######################################################################

install(FILES ${CMAKE_BINARY_DIR}/framework/dmrg/version.h
        DESTINATION include/framework/dmrg COMPONENT headers)

#TODO: improve using GLOB
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

install(DIRECTORY lib/python DESTINATION lib  COMPONENT python
        USE_SOURCE_PERMISSIONS
        FILES_MATCHING PATTERN "*.py" PATTERN "*.pyc"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
