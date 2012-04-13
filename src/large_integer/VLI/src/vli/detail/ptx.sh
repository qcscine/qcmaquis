#!/bin/bash

nvcc -O2  -I ~/VLI_ASM/src/ -I /apps/todi/boost/1.48/gnu_461/include -I ~/VLI_ASM/ -ptx --ptxas-options=-v kernels_gpu.cu 
g++  -E -P -I /apps/eiger/boost_1_46_1/include/ kernels_gpu_asm.hpp | sed  "s/;/;\\`echo -e '\n\r'`/g" > out
