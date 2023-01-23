add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# *** Targets
add_executable(tcDMRG tcDMRG.cpp)
target_link_libraries(tcDMRG ${DMRG_APP_LIBRARIES})

# *** Install
install(TARGETS tcDMRG RUNTIME DESTINATION bin COMPONENT applications)
