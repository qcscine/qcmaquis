add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

set(DMRG_APP_LIBRARIES dmrg_utils ${DMRG_LIBRARIES})

add_executable(mps_overlap mps_overlap.cpp)
target_link_libraries(mps_overlap ${DMRG_APP_LIBRARIES})

add_executable(mps_stats mps_stats.cpp)
target_link_libraries(mps_stats ${DMRG_APP_LIBRARIES})
#set_target_properties(mps_stats PROPERTIES COMPILE_DEFINITIONS "USE_TWOU1")
