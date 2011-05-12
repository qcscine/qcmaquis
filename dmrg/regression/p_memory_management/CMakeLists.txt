add_executable (p_memory_management.test p_memory_management_unit_tests.cpp) 
target_link_libraries (p_memory_management.test ambient boost boost_unit_test_framework alps mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64 pthread) 
add_test(unit_test mpiexec -np 2 ./p_memory_management.test)
