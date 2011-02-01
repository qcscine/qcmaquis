cmake_minimum_required (VERSION 2.6) 
project (AMBIENT_DMRG) 

set(CMAKE_CXX_COMPILER "mpic++") 
#set(CMAKE_VERBOSE_MAKEFILE true)

execute_process(COMMAND rm -rf CMakeCache.txt)

add_subdirectory (ambient) 
add_subdirectory (regression/p_dense_matrix)
