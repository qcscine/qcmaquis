add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

add_executable(mps_overlap mps_overlap.cpp)
target_link_libraries(mps_overlap ${DMRG_APP_LIBRARIES})
