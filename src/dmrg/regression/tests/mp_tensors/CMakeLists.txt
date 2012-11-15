add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

add_executable(observables.test observables.cpp)
target_link_libraries(observables.test ${DMRG_APP_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(observables observables.test)


add_executable(mpsjoin.test mpsjoin.cpp)
target_link_libraries(mpsjoin.test ${DMRG_APP_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(mpsjoin mpsjoin.test)


