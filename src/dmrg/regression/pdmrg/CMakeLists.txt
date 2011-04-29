add_executable (pdmrg.exec ../dmrg/b_u1_pdmrg.cpp ${DMRG_SOURCE_DIR}/mpos/alps_adjacency.cpp)
target_link_libraries (pdmrg.exec ambient  boost boost_unit_test_framework  alps mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64 pthread)
add_test(unit_test mpiexec -np 2 ./pdmrg.exec)
