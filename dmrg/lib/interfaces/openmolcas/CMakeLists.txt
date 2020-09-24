if(NOT CMAKE_Fortran_MODULE_DIRECTORY)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/share/qcmaquis/fortran_mod)
endif()

include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})


message(STATUS "Compiling OpenMolcas interface modules enabled")

set(OPENMOLCAS_INTERFACE_SOURCES
    qcmaquis_interface_cfg.f90
    qcmaquis_interface_utility_routines.f90
    qcmaquis_interface.f90
    qcmaquis_interface_mpssi.f90
    )

set(HDF5_INTERFACE_SOURCES
    hdf5_utils.f90 hdf5_wrapper.c
   )

add_library(qcmaquis-driver ${OPENMOLCAS_INTERFACE_SOURCES})
add_library(qcmaquis-hdf5-interface ${HDF5_INTERFACE_SOURCES})

target_link_libraries(qcmaquis-hdf5-interface ${HDF5_LIBRARIES})

set_target_properties(qcmaquis-driver PROPERTIES COMPILE_FLAGS "-cpp")

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU"
       AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
  target_compile_options(qcmaquis-driver PRIVATE "-fallow-argument-mismatch")
endif()
if (LAPACK_64_BIT)
  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(qcmaquis-driver PROPERTIES COMPILE_FLAGS "-fdefault-integer-8")
    set_source_files_properties(hdf5_utils.f90 PROPERTIES COMPILE_FLAGS "-fdefault-integer-8")
  endif()

  if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set_target_properties(qcmaquis-driver PROPERTIES COMPILE_FLAGS "-i8")
    set_source_files_properties(hdf5_utils.f90 PROPERTIES COMPILE_FLAGS "-i8")
  endif()
endif()

if (BUILD_OPENMOLCAS_MPI)
  target_compile_definitions(qcmaquis-driver PRIVATE "_MOLCAS_MPP_")
  if (GA_INCLUDE_DIR)
    set_target_properties(qcmaquis-driver INCLUDE_DIRECTORIES "${GA_INCLUDE_DIR}")
  endif()
endif()

set_target_properties(qcmaquis-driver PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
set_target_properties(qcmaquis-hdf5-interface PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

install(TARGETS qcmaquis-driver EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
install(TARGETS qcmaquis-hdf5-interface EXPORT DMRGTargets COMPONENT libraries DESTINATION lib)
install(PROGRAMS tools/qcm_checkpoint_rename.py DESTINATION bin COMPONENT tools)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION share/qcmaquis)

message(STATUS "${PROJECT_BINARY_DIR}")
export(TARGETS qcmaquis-driver APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
export(TARGETS qcmaquis-hdf5-interface APPEND FILE "${PROJECT_BINARY_DIR}/DMRGTargets.cmake")
