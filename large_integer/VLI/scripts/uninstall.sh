#!/bin/bash

#global variables
PATH_SRC=${PWD}/..

delete(){
    cd ${PATH_SRC}
    rm -rf build_*
}

echo "start cleaning"
delete
echo "end cleaning"
