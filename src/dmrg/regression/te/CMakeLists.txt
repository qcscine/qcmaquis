
set(DMRG_APP_LIBRARIES ${DMRG_LIBRARIES})

add_executable (block_to_mpo.test ${CMAKE_SOURCE_DIR}/utils/utils.cpp block_to_mpo.cpp)
target_link_libraries(block_to_mpo.test ${DMRG_APP_LIBRARIES})
add_test(block_to_mpo ${CMAKE_CURRENT_BINARY_DIR}/block_to_mpo.test)


add_executable (make_mpo_exp.test ${CMAKE_SOURCE_DIR}/utils/utils.cpp make_mpo_exp.cpp)
target_link_libraries(make_mpo_exp.test ${DMRG_APP_LIBRARIES})
add_test(make_mpo_exp ${CMAKE_CURRENT_BINARY_DIR}/make_mpo_exp.test)
