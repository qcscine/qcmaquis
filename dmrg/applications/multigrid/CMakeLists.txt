
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)


# *** Targets
add_executable(multigrid_dmrg_u1 multigrid.cpp)
target_link_libraries(multigrid_dmrg_u1 ${DMRG_APP_LIBRARIES})
add_executable(multigrid_dmrg_2u1 multigrid.cpp)
target_link_libraries(multigrid_dmrg_2u1 ${DMRG_APP_LIBRARIES})

add_executable(multigrid_measure_u1 multigrid.cpp)
target_link_libraries(multigrid_measure_u1 ${DMRG_APP_LIBRARIES})
add_executable(multigrid_measure_2u1 multigrid.cpp)
target_link_libraries(multigrid_measure_2u1 ${DMRG_APP_LIBRARIES})

set_target_properties(multigrid_measure_u1 multigrid_measure_2u1 PROPERTIES COMPILE_DEFINITIONS MEASURE_ONLY)
set_target_properties(multigrid_dmrg_2u1 multigrid_measure_2u1 PROPERTIES COMPILE_DEFINITIONS UseTwoU1)



# *** Install
install(TARGETS multigrid_dmrg_u1 multigrid_dmrg_2u1 multigrid_measure_u1 multigrid_measure_2u1 RUNTIME DESTINATION bin COMPONENT applications)
