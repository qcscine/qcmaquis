add_executable (dmrg.exec b_u1_dmrg.cpp ${DMRG_SOURCE_DIR}/mpos/alps_adjacency.cpp)
target_link_libraries (dmrg.exec boost alps mkl_intel_lp64 mkl_sequential mkl_core pthread )
