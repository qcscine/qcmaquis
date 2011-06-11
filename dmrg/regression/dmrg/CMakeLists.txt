add_definitions(-DHAVE_ALPS_HDF5)
add_executable (dmrg.exec b_u1_dmrg.cpp ${DMRG_SOURCE_DIR}/mpos/alps_adjacency.cpp)
target_link_libraries (dmrg.exec ${DMRG_LIBRARIES})
