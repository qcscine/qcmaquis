#!/bin/bash -l

PATH_SRC=${PWD}/..
MACHINE_CONFIG="castor"

run_dashboad(){
    cd ${PATH_SRC}/build_${MACHINE_CONFIG}_${COMPILER_NAME}
    make clean
    make
    ctest -S Dashboards/site.cmake
}

module load cmake
COMPILER_NAME="icc"
module load intel
run_dashboad

COMPILER_NAME="gcc"
module load gcc
run_dashboad
