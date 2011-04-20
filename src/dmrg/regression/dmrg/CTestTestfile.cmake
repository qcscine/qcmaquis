# CMake generated Testfile for 
# Source directory: /users/ewartt/DMRG/regression/dmrg
# Build directory: /users/ewartt/DMRG/regression/dmrg
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(unit_test "mpiexec" "-np" "2" "./p_dense_matrix.test")
