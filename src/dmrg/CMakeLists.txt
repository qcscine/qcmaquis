cmake_minimum_required (VERSION 2.6) 
project (DMRG) 
enable_testing()
set (DMRG_VERSION_MAJOR 0)
set (DMRG_VERSION_MINOR 1) 
# use for make debugging
# set(CMAKE_VERBOSE_MAKEFILE true)
macro(add arg) 
if(NOT DEFINED ${arg})
add_subdirectory (${arg}) 
set( ${arg} true)
endif(NOT DEFINED ${arg}) 
endmacro(add)
macro(target arg) 
set( ${arg} true) 
endmacro(target)

# possible build targets:
target(ambient_dmrg)

# description of targets:
if(ambient_dmrg)
# note: you can regenerate operation.pp.h/hpp by invoking smth like the following:
# g++ ambient/core/operation/operation.pp.h.template.h   -P -E -Iboost -I. > operation.pp.h
# g++ ambient/core/operation/operation.pp.hpp.template.h -P -E -Iboost -I. > operation.pp.hpp
# g++ ambient/interface/push.pp.hpp.template.h           -P -E -Iboost -I. > push.pp.hpp
    add_definitions(-O0 -g -w)
    add_definitions(-c -DSCALAPACK -DBOOST_TEST_DYN_LINK -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DCPU_FREQ=2.5e9)
    include_directories(${DMRG_SOURCE_DIR}
                        /apps/eiger/Intel-MPI-4.0/intel64/include
                        /users/bauerb/ALPS-eiger/include
                        /users/bauerb/alps_svn/boost_1_45_0
                        /apps/eiger/Intel-FOR-11.1/mkl/include)
    link_directories   (${DMRG_BINARY_DIR}/ambient
                        /users/bauerb/ALPS-eiger/lib
                        /apps/eiger/Intel-MPI-4.0/intel64/lib
                        /apps/eiger/Intel-FOR-11.1/mkl/lib/em64t)
    set(CMAKE_CXX_COMPILER "g++") 
    #add(regression/dmrg)

    add_definitions(-DMPI_PARALLEL)
    set(CMAKE_CXX_COMPILER "mpiicpc") 

    add(ambient) 
    add(regression/p_dense_matrix)
    add(regression/p_diagonal_matrix)
    add(regression/p_reshape)
    add(regression/p_dense_matrix_benchmark)
    add(regression/pdmrg)
endif(ambient_dmrg)


