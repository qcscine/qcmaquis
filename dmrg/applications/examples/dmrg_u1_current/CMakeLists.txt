add_definitions(-DHAVE_ALPS_HDF5 -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)
enable_omp_if_found()
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

set(DMRG_APP_LIBRARIES dmrg_utils dmrg_models ${DMRG_LIBRARIES})

# *** Targets
add_executable(dmrg_u1_current dmrg_u1_current.cpp)
target_link_libraries(dmrg_u1_current ${DMRG_APP_LIBRARIES})

