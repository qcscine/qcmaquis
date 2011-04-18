
add_definitions(-c -O0 -g -DSCALAPACK -DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DCPU_FREQ=2.5e9 -DHAVE_ALPS_HDF5 -DMPI_PARALLEL)

include_directories (${DMRG_SOURCE_DIR}
                     /apps/eiger/Intel-MPI-4.0/intel64/include
                     /apps/eiger/boost-1.43/include
                     /users/bauerb/ALPS-eiger/include
                     /apps/eiger/Intel-FOR-11.1/mkl/include
                     /apps/eiger/hdf5/1.8.7/gcc-ifort/include)

link_directories (${DMRG_BINARY_DIR}/ambient
                  /apps/eiger/hdf5/1.8.7/gcc-ifort/lib
                  /cluster/home/phys/bauerb/ALPS/lib
                  /apps/eiger/boost-1.43/lib
                  /users/bauerb/ALPS-eiger/lib
                  /apps/eiger/Intel-MPI-4.0/intel64/lib
                  /apps/eiger/Intel-FOR-11.1/mkl/lib/em64t)

add_executable (p_dense_matrix.benchmark p_dense_matrix_benchmark.cpp) 
target_link_libraries (p_dense_matrix.benchmark ambient boost boost_system boost_unit_test_framework alps hdf5  mkl_scalapack_lp64 mkl_solver_lp64_sequential mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64 pthread ) 


add_test(unit_test mpiexec -np 2 ./p_dense_matrix.benchmark)
