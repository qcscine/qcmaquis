
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})
if(USE_MTM)
  list(INSERT DMRG_APP_LIBRARIES 0 mtm)
endif(USE_MTM)


# *** Targets
add_executable(te_nn_u1 te_nn.cpp)
target_link_libraries(te_nn_u1 ${DMRG_APP_LIBRARIES})

add_executable(te_nn_2u1 te_nn.cpp)
target_link_libraries(te_nn_2u1 ${DMRG_APP_LIBRARIES})

add_executable(te_u1 te.cpp)
target_link_libraries(te_u1 ${DMRG_APP_LIBRARIES})

add_executable(te_2u1 te.cpp)
target_link_libraries(te_2u1 ${DMRG_APP_LIBRARIES})


set_target_properties(te_nn_2u1 te_2u1 PROPERTIES COMPILE_DEFINITIONS UseTwoU1)



install(TARGETS te_nn_u1 te_nn_2u1 te_u1 te_2u1 RUNTIME DESTINATION bin COMPONENT applications)

