add_executable (p_diagonal_matrix.test p_diagonal_matrix_unit_tests.cpp) 
target_link_libraries (p_diagonal_matrix.test ambient boost boost_unit_test_framework alps mkl_scalapack_lp64 mkl_solver_lp64_sequential mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64 pthread) 
add_test(unit_test mpiexec -np 2 ./p_diagonal_matrix.test)
