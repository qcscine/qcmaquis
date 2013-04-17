compile_with_ambient()
add_executable (pdmrg.exec ../dmrg/b_u1_pdmrg.cpp
                           ${CMAKE_SOURCE_DIR}/framework/dmrg/utils/utils.cpp
                           ${CMAKE_SOURCE_DIR}/framework/dmrg/deprecated/mpos/alps_adjacency.cpp)
target_link_libraries (pdmrg.exec alps ${AMBIENT_LIBRARIES} boost ${BLAS_LAPACK_LIBRARIES})
add_test(unit_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_POSTFLAGS} ${MPIEXEC_PREFLAGS} ./pdmrg.exec)

