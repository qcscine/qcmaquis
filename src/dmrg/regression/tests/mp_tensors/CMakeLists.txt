add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

add_executable(observables.test observables.cpp)
target_link_libraries(observables.test ${DMRG_APP_LIBRARIES})

add_test(observables observables.test)


add_executable(mpsjoin.test mpsjoin.cpp)
target_link_libraries(mpsjoin.test ${DMRG_APP_LIBRARIES})

add_test(mpsjoin mpsjoin.test)


add_executable(coherent_init.test coherent_init.cpp)
target_link_libraries(coherent_init.test ${DMRG_APP_LIBRARIES})

add_test(coherent_init coherent_init.test)


add_executable(apply_op.test apply_op.cpp)
target_link_libraries(apply_op.test dmrg_models ${DMRG_APP_LIBRARIES})

add_test(apply_op apply_op.test)


add_executable(super_mpo.test super_mpo.cpp)
target_link_libraries(super_mpo.test ${DMRG_APP_LIBRARIES})

add_test(super_mpo super_mpo.test)


