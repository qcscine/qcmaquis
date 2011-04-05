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
    set(CMAKE_CXX_COMPILER "mpic++") 
    add_definitions(-O3)
    add(ambient) 
    #add(regression/p_dense_matrix)
    add(regression/p_diagonal_matrix)
endif(ambient_dmrg)


