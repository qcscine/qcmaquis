
set(DMRG_APP_LIBRARIES ${DMRG_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR}/applications/time_evolution)

add_executable (block_to_mpo.test ${DMRG_FRAMEWORK_DIR}/dmrg/utils/utils.cpp block_to_mpo.cpp)
target_link_libraries(block_to_mpo.test ${DMRG_APP_LIBRARIES})
add_test(NAME block_to_mpo
         # CONFIGURATIONS Debug    
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/block_to_mpo.test
         )


add_executable (make_mpo_exp.test ${DMRG_FRAMEWORK_DIR}/dmrg/utils/utils.cpp make_mpo_exp.cpp)
target_link_libraries(make_mpo_exp.test ${DMRG_APP_LIBRARIES})
add_test(NAME make_mpo_exp
         # CONFIGURATIONS Debug    
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/make_mpo_exp.test
         )
