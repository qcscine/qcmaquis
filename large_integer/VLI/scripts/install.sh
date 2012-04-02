#!/bin/bash -l

#paths 
PATH_SRC=${PWD}/..
PATH_SCRIPT=${PWD}
export BOOST_ROOT=/apps/eiger/boost_1_46_1/
#some export compiler and boost
create_directory(){
    cd ${PATH_SRC}
    mkdir build_${MACHINE_CONFIG}_${COMPILER_NAME}
}

create_vli(){
    cd ${PATH_SRC}/build_${MACHINE_CONFIG}_${COMPILER_NAME}
    create_dashboards
    echo "Machine config \" ${PATH_SCRIPT}/${MACHINE_CONFIG}_${COMPILER_NAME}.cmake\" "
    cmake -DMACHINE_CONFIG=${PATH_SCRIPT}/${MACHINE_CONFIG}_${COMPILER_NAME}.cmake -DVLI_TESTS=ON -DVLI_MAIN=OFF ..
}

create_dashboards(){
    mkdir Dashboards
    echo "set(PREDEFINED_CTEST_SITE \"`uname -n`\")"                    >  ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BUILD_NAME \"${COMPILER_NAME}\")"        >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_SOURCE_DIRECTORY \"${PATH_SRC}\")"       >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BINARY_DIRECTORY \"${PATH_SRC}/build_${MACHINE_CONFIG}_${COMPILER_NAME}\")" >> ./Dashboards/site.cmake
    cat ../Dashboards/site.cmake >> ./Dashboards/site.cmake
    cp ../Dashboards/cmake_common.cmake ./Dashboards/
}

create_distrib(){
    echo "Start creation VLI"
    create_directory
    create_vli
    echo "End creation VLI"
}

module load cmake
MACHINE_CONFIG="castor"
COMPILER_NAME="icc"
module load intel
create_distrib

COMPILER_NAME="gcc"
module load gcc
create_distrib


