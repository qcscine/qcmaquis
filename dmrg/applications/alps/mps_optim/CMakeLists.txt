
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

# *** Targets
add_executable(mps_optim run_dmrg.cpp scheduler.cpp mps_optim.cpp)
target_link_libraries(mps_optim ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS mps_optim RUNTIME DESTINATION bin COMPONENT applications)
