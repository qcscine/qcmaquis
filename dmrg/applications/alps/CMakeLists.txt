
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# *** Libraries

add_library(mps_pscan STATIC libpscan/scheduler.cpp libpscan/options.cpp)



set(MYAPP_LIBRARIES mps_pscan dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# *** Targets
add_executable(mps_optim mps_optim/main.cpp mps_optim/run_sim.cpp)
target_link_libraries(mps_optim ${MYAPP_LIBRARIES})

add_executable(mps_evolve mps_evolve/main.cpp mps_evolve/run_sim.cpp)
target_link_libraries(mps_evolve ${MYAPP_LIBRARIES})

# *** Install
install(TARGETS mps_optim RUNTIME DESTINATION bin COMPONENT applications)
