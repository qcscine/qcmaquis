

add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)

if(USE_AMBIENT)
  compile_with_ambient()
  list(INSERT DMRG_APP_LIBRARIES 0 ambient ${MPI_CXX_LIBRARIES})
endif(USE_AMBIENT)  

# *** Configure symmetries
configure_symm_file("mg_dmrg_symm/mg_dmrg_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/mg_dmrg_symm/mg_dmrg_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND DMRG_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/mg_dmrg_symm/mg_dmrg_{SYMM}.cpp")

configure_symm_file("mg_meas_symm/mg_meas_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/mg_meas_symm/mg_meas_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND MEASURE_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/mg_meas_symm/mg_meas_{SYMM}.cpp")


# *** Targets
add_executable(mg_dmrg ${DMRG_SYMM_SOURCES} mg_dmrg.cpp)
target_link_libraries(mg_dmrg ${DMRG_APP_LIBRARIES})

add_executable(mg_meas ${MEASURE_SYMM_SOURCES} mg_meas.cpp)
target_link_libraries(mg_meas ${DMRG_APP_LIBRARIES})


# *** Install
install(TARGETS mg_dmrg mg_meas RUNTIME DESTINATION bin COMPONENT applications)


if(BUILD_OLD_TARGETS)
    # *** Targets
    add_executable(multigrid_dmrg_u1 multigrid.cpp)
    target_link_libraries(multigrid_dmrg_u1 ${DMRG_APP_LIBRARIES})
    add_executable(multigrid_dmrg_2u1 multigrid.cpp)
    target_link_libraries(multigrid_dmrg_2u1 ${DMRG_APP_LIBRARIES})

    add_executable(multigrid_measure_u1 multigrid.cpp)
    target_link_libraries(multigrid_measure_u1 ${DMRG_APP_LIBRARIES})
    add_executable(multigrid_measure_2u1 multigrid.cpp)
    target_link_libraries(multigrid_measure_2u1 ${DMRG_APP_LIBRARIES})

    set_property(TARGET multigrid_measure_u1 multigrid_measure_2u1 APPEND PROPERTY COMPILE_DEFINITIONS MEASURE_ONLY)
    set_property(TARGET multigrid_dmrg_2u1 multigrid_measure_2u1 APPEND PROPERTY COMPILE_DEFINITIONS UseTwoU1)

    # *** Install
    install(TARGETS multigrid_dmrg_u1 multigrid_dmrg_2u1 multigrid_measure_u1 multigrid_measure_2u1 RUNTIME DESTINATION bin COMPONENT applications)
endif(BUILD_OLD_TARGETS)
