
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES maquis_dmrg-shared dmrg_models dmrg_utils ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)

# *** Targets
add_executable(dmrg_meas measure.cpp)
target_link_libraries(dmrg_meas ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS dmrg_meas RUNTIME DESTINATION bin COMPONENT applications)
