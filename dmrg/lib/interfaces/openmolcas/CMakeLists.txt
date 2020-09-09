project (MAQUIS_DMRG CXX C Fortran)

if(NOT CMAKE_Fortran_MODULE_DIRECTORY)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/mod)
endif()

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})


message(STATUS "Compiling OpenMolcas interface modules")
message(STATUS "-- Module directory: ${CMAKE_Fortran_MODULE_DIRECTORY}" )

set(OPENMOLCAS_INTERFACE_SOURCES
    qcmaquis_interface_cfg.f90
    qcmaquis_interface_utility_routines.f90
    qcmaquis_interface.f90
    qcmaquis_interface_mpssi.f90
    )

set(HDF5_INTERFACE_SOURCES
    hdf5_utils.f90 hdf5_wrapper.c
   )

add_library(qcmaquis-driver-utils ${OPENMOLCAS_INTERFACE_SOURCES})
add_library(HDF5interface ${HDF5_INTERFACE_SOURCES})

target_link_libraries(HDF5interface ${HDF5_LIBRARIES})

set_target_properties(qcmaquis-driver-utils PROPERTIES COMPILE_FLAGS "-cpp")

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU"
       AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
  target_compile_options(qcmaquis-driver-utils PRIVATE "-fallow-argument-mismatch")
endif()
if (LAPACK_64_BIT)
  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(qcmaquis-driver-utils PROPERTIES COMPILE_FLAGS "-fdefault-integer-8")
    set_source_files_properties(hdf5_utils.f90 PROPERTIES COMPILE_FLAGS "-fdefault-integer-8")
  endif()

  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set_target_properties(qcmaquis-driver-utils PROPERTIES COMPILE_FLAGS "-i8")
    set_source_files_properties(hdf5_utils.f90 PROPERTIES COMPILE_FLAGS "-i8")
  endif()
endif()

if (BUILD_OPENMOLCAS_MPI)
  target_compile_definitions(qcmaquis-driver-utils PRIVATE "_MOLCAS_MPP_")
  if (GA_INCLUDE_DIR)
    set_target_properties(qcmaquis-driver-utils INCLUDE_DIRECTORIES "${GA_INCLUDE_DIR}")
  endif()
endif()

set_target_properties(qcmaquis-driver-utils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
set_target_properties(HDF5interface PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

install(TARGETS qcmaquis-driver-utils EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
install(TARGETS HDF5interface EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)

export(TARGETS qcmaquis-driver-utils APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
export(TARGETS HDF5interface APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
