cmake_minimum_required(VERSION 3.16.3)
project(MAQUIS_DMRG C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(QCMAQUIS_TESTS)
  message(STATUS "Enabling tests")
  enable_testing()
endif()

add_subdirectory(src)
