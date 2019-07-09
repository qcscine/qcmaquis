add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

add_executable(custom_model.test custom_model.cpp)
target_link_libraries(custom_model.test ${DMRG_APP_LIBRARIES})

add_test(custom_model custom_model.test)
