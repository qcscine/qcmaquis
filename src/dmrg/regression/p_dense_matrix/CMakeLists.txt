
add_definitions(-c -O3 -DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DCPU_FREQ=2.5e9 -DHAVE_ALPS_HDF5 -DMPI_PARALLEL)

include_directories (${DMRG_SOURCE_DIR}
                     /Users/akosenko/workspace/build/boost_1_45_0
                     /Users/akosenko/workspace/build/alps_14.12.2010/include
                     /Users/akosenko/workspace/build/hdf5_1.8.5/include)

link_directories (${DMRG_BINARY_DIR}/ambient
                  /Users/akosenko/workspace/build/hdf5_1.8.5/lib
                  /Users/akosenko/workspace/build/alps_14.12.2010/lib
                  /opt/local/lib
                  /opt/boost-1.44-binary/lib) 

add_executable (p_dense_matrix.test p_dense_matrix_unit_tests.cpp) 
target_link_libraries (p_dense_matrix.test ambient boost alps hdf5 arpack z sz boost_unit_test_framework) 

add_test(unit_test mpiexec -np 2 ./p_dense_matrix.test)
