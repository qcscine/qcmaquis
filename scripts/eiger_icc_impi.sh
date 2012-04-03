#!/bin/bash

SITE=eiger # the name of the site on the dashboard # or `uname -n`
ROOT_DIR=~/maquis2012/src
ALPS_ROOT_DIR=/project/h07/ALPS_INTEL_EIGER
BOOST_BINDINGS_INCLUDE=/project/h07/ALPS_INTEL_EIGER/include
MACHINE_CMAKE_CONFIG='eiger_intel.cmake'

COMPILER=iccxe
COMPILER_VERSION=2011
MPI_WRAPPER=impi
MPI_WRAPPER_VERSION=4.0
BOOST=boost
BOOST_VERSION=1.46.1

export CXX=icpc
export CC=icc
export BOOST_ROOT=/apps/eiger/boost_1_46_1

source common.sh 

module load ${COMPILER}/${COMPILER_VERSION} ${BOOST}/${BOOST_VERSION} ${MPI_WRAPPER}/${MPI_WRAPPER_VERSION} 
create_build_tree
