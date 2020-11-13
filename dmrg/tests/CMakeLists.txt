add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# *** Targets
add_executable(test1 test1.cpp)
target_link_libraries(test1 ${DMRG_APP_LIBRARIES})

add_executable(test2 test2.cpp)
target_link_libraries(test2 ${DMRG_APP_LIBRARIES})

add_executable(test_integral_map test_integral_map.cpp)
target_link_libraries(test_integral_map ${DMRG_APP_LIBRARIES})

add_executable(test_rel test_rel.cpp)
target_link_libraries(test_rel ${DMRG_APP_LIBRARIES})

add_executable(test_4rdm_save test_4rdm_save.cpp)
target_link_libraries(test_4rdm_save ${DMRG_APP_LIBRARIES})

add_executable(test_mps_rotate test_mps_rotate.cpp)
target_link_libraries(test_mps_rotate ${DMRG_APP_LIBRARIES})

add_executable(test_wigner test_wigner.cpp)
target_link_libraries(test_wigner ${DMRG_APP_LIBRARIES})

# -- Unit tests --
add_executable(test_block_matrix block_matrix/block_matrix.cpp)
target_link_libraries(test_block_matrix ${DMRG_APP_LIBRARIES})
