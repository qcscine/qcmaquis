configure_symm_file("simulation_symm/simulation_tpl.cpp.in"
                    "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
set(MPS_EIGENSTATE_RUN_SOURCES)
get_symm_files(APPEND MPS_EIGENSTATE_RUN_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")
# TODO: change to add_library( OBJECT ) once CMake 2.8.8 is the minimum required
add_library(mps_eigenstate_run STATIC ${MPS_EIGENSTATE_RUN_SOURCES})


add_executable(mps_optim mps_optim.cpp run_optim.cpp)
target_link_libraries(mps_optim mps_eigenstate_run ${MYAPP_LIBRARIES})

add_executable(mps_meas mps_meas.cpp run_measure.cpp)
target_link_libraries(mps_meas mps_eigenstate_run ${MYAPP_LIBRARIES})


install(TARGETS mps_optim mps_meas RUNTIME DESTINATION bin COMPONENT applications)
