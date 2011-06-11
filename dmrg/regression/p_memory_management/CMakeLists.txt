enable_ambient()
add_executable (p_memory_management.test p_memory_management_unit_tests.cpp) 
target_link_libraries (p_memory_management.test ${AMBIENT_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} boost_unit_test_framework) 
add_test(unit_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_POSTFLAGS} ${MPIEXEC_PREFLAGS} ./p_memory_management.test)
