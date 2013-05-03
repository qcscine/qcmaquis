add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

add_executable(super_meas.test super_meas.cpp)
target_link_libraries(super_meas.test ${DMRG_APP_LIBRARIES})

add_test(super_meas super_meas.test)

