add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

add_executable(qrdecomp.test qrdecomp.cpp)
target_link_libraries(qrdecomp.test ${DMRG_APP_LIBRARIES} ${MAQUIS_Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(qrdecomp qrdecomp.test)
