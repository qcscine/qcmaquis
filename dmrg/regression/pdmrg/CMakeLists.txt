set(CMAKE_CXX_COMPILER "mpic++") 
add_definitions(-c -w0  -O0 -g  -DSCALAPACK -DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DCPU_FREQ=2.2e9 -DHAVE_ALPS_HDF5 -DMPI_PARALLEL)

include_directories (${DMRG_SOURCE_DIR}
	  		/users/bauerb/ALPS-eiger/include 
                        /users/bauerb/alps_svn/boost_1_45_0	  		
                        /apps/eiger/hdf5/1.8.7/gcc-ifort/include
                        /apps/eiger/Intel-FOR-11.1/mkl/include)

link_directories (${DMRG_BINARY_DIR}/ambient
                        /apps/eiger/hdf5/1.8.7/gcc-ifort/lib
                        /users/bauerb/ALPS-eiger/lib  
                        /apps/eiger/Intel-FOR-11.1/mkl/lib/em64t)

add_executable (dmrg_mpi.exec ../dmrg/b_u1_dmrg.cpp ${DMRG_SOURCE_DIR}/mpos/alps_adjacency.cpp)

#target_link_libraries (dmrg_mpi.exec  boost alps mkl_intel_lp64 mkl_sequential mkl_core pthread)
target_link_libraries (dmrg_mpi.exec ambient  boost boost_unit_test_framework  alps mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_openmpi_lp64 pthread )

add_test(unit_test mpiexec -np 2 ./p_dense_matrix.test)
