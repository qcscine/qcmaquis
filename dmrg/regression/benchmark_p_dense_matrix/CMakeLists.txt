
add_definitions(-c -O0 -g -DSCALAPACK -DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DCPU_FREQ=2.5e9 -DHAVE_ALPS_HDF5 -DMPI_PARALLEL)

include_directories (${DMRG_SOURCE_DIR}
/cluster/home/phys/ewartt/BOOST/boost_1_45_0/build/include
                     /cluster/home/phys/bauerb/ALPS/include
                     /cluster/apps/intel/mkl/10.1.1.019/include)

link_directories (${DMRG_BINARY_DIR}/ambient
/cluster/home/phys/ewartt/BOOST/boost_1_45_0/build/lib
                  /cluster/home/phys/bauerb/ALPS/lib
                  /cluster/apps/intel/mkl/10.1.1.019/lib)

add_executable (bench.exe benchmark_p_dense_matrix.cpp) 
target_link_libraries (bench.exe ambient boost alps hdf5 boost_unit_test_framework  mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_openmpi_lp64 pthread ) 

add_test(unit_test mpiexec -np 2 ./p_dense_matrix.test)
