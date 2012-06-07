#!/bin/bash

g++  -E -P -I /usr/local/include/boost/include kernels_gpu_asm.hpp | sed  "s/;/;\\`echo -e '\n\r'`/g" > out
