add_executable (p_dense_matrix.benchmark p_dense_matrix_benchmark.cpp) 
target_link_libraries (p_dense_matrix.benchmark ambient boost boost_unit_test_framework alps mkl_scalapack_lp64 mkl_solver_lp64_sequential mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64 pthread) 
add_test(unit_test mpiexec -np 2 ./p_dense_matrix.benchmark)
