#!/bin/bash

## targets ##

AMBIENT=ambient
DMRG=dmrg
TYPES=types
BUILD_DIR=${SITE}_${COMPILER}_${MPI_WRAPPER}

## functions ##

print_configure_begin(){
    echo " ----------------------- "
    echo " configure build tree    "
    echo " ${BUILD_DIR}            "
    echo " ----------------------- "
}

print_configure_end(){
    echo " ----------------------- "
    echo " build configure is done "
    echo " ${BUILD_DIR}            "
    echo " ----------------------- "
}

add_ambient(){
    mkdir ${ROOT_DIR}/${AMBIENT}/${BUILD_DIR}
    cd ${ROOT_DIR}/${AMBIENT}/${BUILD_DIR}
    cmake -DMACHINE_CONFIG=${ROOT_DIR}/${DMRG}/config_machines/${MACHINE_CMAKE_CONFIG} ..
}

remove_ambient(){
    rm -rf ${ROOT_DIR}/${AMBIENT}/${BUILD_DIR}
}

add_dmrg(){
    mkdir ${ROOT_DIR}/${DMRG}/${BUILD_DIR}
    cd ${ROOT_DIR}/${DMRG}/${BUILD_DIR}
    cmake -DMACHINE_CONFIG=${ROOT_DIR}/${DMRG}/config_machines/${MACHINE_CMAKE_CONFIG} -DALPS_ROOT_DIR=${ALPS_ROOT_DIR} -DBUILD_REGRESSION=ON -DBUILD_AMBIENT=ON -DUSE_AMBIENT=ON .. 
}

remove_dmrg(){
    rm -rf ${ROOT_DIR}/${DMRG}/${BUILD_DIR}
}

add_types(){
    mkdir ${ROOT_DIR}/${TYPES}/${BUILD_DIR}
    cd ${ROOT_DIR}/${TYPES}/${BUILD_DIR}
    use_dashboards
    cmake -DMACHINE_CONFIG=${ROOT_DIR}/${DMRG}/config_machines/${MACHINE_CMAKE_CONFIG} -DBOOST_BINDINGS_INCLUDE=${BOOST_BINDINGS_INCLUDE} -DENABLE_PARALLEL=ON -DBUILD_AMBIENT=ON -DENABLE_REGRESSION_FUNCTIONAL=ON ..
}

remove_types(){
    rm -rf ${ROOT_DIR}/${TYPES}/${BUILD_DIR}
}

use_dashboards(){
    mkdir Dashboards
    echo "set(PREDEFINED_CTEST_SITE \"${SITE}\")"                                       >  ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BUILD_NAME \"${COMPILER}_${MPI_WRAPPER}\")"              >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_SOURCE_DIRECTORY \"${ROOT_DIR}/${TYPES}\")"              >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BINARY_DIRECTORY \"${ROOT_DIR}/${TYPES}/${BUILD_DIR}\")" >> ./Dashboards/site.cmake
    cat ../Dashboards/site.cmake                                                        >> ./Dashboards/site.cmake
    cp ../Dashboards/cmake_common.cmake ./Dashboards/
}

create_build_tree(){
    print_configure_begin
    add_ambient
    add_dmrg
    add_types
    print_configure_end
}

delete_buid_tree(){
    remove_ambient
    remove_dmrg
    remove_types
}
