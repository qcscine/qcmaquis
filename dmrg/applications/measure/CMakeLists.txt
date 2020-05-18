
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_models dmrg_utils mpi_interface ${DMRG_LIBRARIES})

if(USE_AMBIENT)
  compile_with_ambient()
endif(USE_AMBIENT)  

# *** Configure symmetries
configure_symm_file("simulation_symm/simulation_tpl.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp" CMAKE_SYMM_GROUP)
get_symm_files(APPEND MEASURE_SYMM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/simulation_symm/simulation_{SYMM}.cpp")


# *** Targets
add_executable(dmrg_meas ${MEASURE_SYMM_SOURCES} measure.cpp)
target_link_libraries(dmrg_meas ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS dmrg_meas RUNTIME DESTINATION bin COMPONENT applications)
