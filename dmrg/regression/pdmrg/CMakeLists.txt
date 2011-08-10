enable_ambient()
add_executable (pdmrg.exec ../dmrg/b_u1_pdmrg.cpp
${DMRG_SOURCE_DIR}/utils/utils.cpp ${DMRG_SOURCE_DIR}/mpos/alps_adjacency.cpp)
target_link_libraries (pdmrg.exec ${AMBIENT_SCALAPACK_LIBRARIES} ${BLAS_LAPACK_LIBRARIES} alps boost)
add_test(unit_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_POSTFLAGS} ${MPIEXEC_PREFLAGS} ./pdmrg.exec)
