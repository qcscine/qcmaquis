add_definitions(-DHAVE_ALPS_HDF5)
add_executable (dmrg.exec b_u1_dmrg.cpp
                          ${CMAKE_SOURCE_DIR}/framework/dmrg/utils/utils.cpp
                          ${CMAKE_SOURCE_DIR}/framework/dmrg/utils/stream_storage.cpp
                          ${CMAKE_SOURCE_DIR}/framework/dmrg/deprecated/mpos/alps_adjacency.cpp)

target_link_libraries (dmrg.exec ${DMRG_LIBRARIES})
