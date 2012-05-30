
set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

add_executable(block_matrix block_matrix.cpp)
target_link_libraries(block_matrix ${DMRG_APP_LIBRARIES})

add_test(block_matrix block_matrix) 
