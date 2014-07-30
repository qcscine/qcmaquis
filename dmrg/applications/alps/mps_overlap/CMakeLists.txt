set(MPS_OPTIM_SOURCES main.cpp run_sim.cpp)

add_executable(mps_overlap main.cpp)
target_link_libraries(mps_overlap ${MYAPP_LIBRARIES})

install(TARGETS mps_overlap RUNTIME DESTINATION bin COMPONENT applications)
