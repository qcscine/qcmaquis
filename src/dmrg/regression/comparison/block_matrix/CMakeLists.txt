
set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

message(${DMRG_APP_LIBRARIES})

add_executable(reshape reshape.cpp)
target_link_libraries(reshape  ${DMRG_APP_LIBRARIES})

add_test(reshape reshape) 
