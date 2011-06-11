enable_ambient()
add_executable (p_dense_matrix.test p_dense_matrix_unit_tests.cpp) 
target_link_libraries (p_dense_matrix.test ${AMBIENT_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} boost_unit_test_framework) 
add_test(unit_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_POSTFLAGS} ${MPIEXEC_PREFLAGS} ./p_dense_matrix.test)
