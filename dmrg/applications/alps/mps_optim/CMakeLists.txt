configure_symm_file("simulation_symm/simulation_tpl.cpp.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
set(MPS_OPTIM_SOURCES main.cpp run_sim.cpp)
get_symm_files(APPEND MPS_OPTIM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")

add_executable(mps_optim ${MPS_OPTIM_SOURCES})
target_link_libraries(mps_optim ${MYAPP_LIBRARIES})

install(TARGETS mps_optim RUNTIME DESTINATION bin COMPONENT applications)
