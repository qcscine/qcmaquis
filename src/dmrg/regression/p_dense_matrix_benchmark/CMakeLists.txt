enable_ambient()
add_executable (p_dense_matrix.benchmark ${PROJECT_SOURCE_DIR}/utils/utils.cpp p_dense_matrix_benchmark.cpp) 
target_link_libraries (p_dense_matrix.benchmark ${AMBIENT_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} boost) 
add_test(unit_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_POSTFLAGS} ${MPIEXEC_PREFLAGS} ./p_dense_matrix.benchmark)
