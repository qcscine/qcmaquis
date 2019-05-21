add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES maquis_dmrg dmrg_models dmrg_utils ${DMRG_LIBRARIES})

# *** Targets
add_executable(test1 test1.cpp)
target_link_libraries(test1 ${DMRG_APP_LIBRARIES})

add_executable(test2 test2.cpp)
target_link_libraries(test2 ${DMRG_APP_LIBRARIES})
