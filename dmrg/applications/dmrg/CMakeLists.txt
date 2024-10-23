
add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# *** Targets
add_executable(dmrg dmrg.cpp)
target_link_libraries(dmrg ${DMRG_APP_LIBRARIES})

add_executable(dmrgIPI dmrg_ipi.cpp)
target_link_libraries(dmrgIPI ${DMRG_APP_LIBRARIES})

add_executable(fww-mpo fww-mpo.cpp)
target_link_libraries(fww-mpo ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS dmrgIPI RUNTIME DESTINATION bin COMPONENT applications)
