#all: int64-version gmp-version vli-version
CC=g++
#CC=llvm-g++

BOOST_PATH=${HOME}/boost_1_45_0
VLI_PATH=${HOME}/vli_install
CUDA_PATH=/cluster/apps/cuda/4.0
GMP_PATH=${HOME}/opt/gmp

CXXFLAGS=-Wall -O3 -funroll-loops -fopenmp
#COMMON_DEFS= -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG 
COMMON_DEFS= -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DVLI_MAX_ORDER=4 -DVECTOR_MAX_ORDER=512 

all: vlicpu gmp

vligpu:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) -I$(BOOST_PATH)/include minimal_high_t_expansion.cpp -DUSE_VLI_INTEGERS_GPU -I$(VLI_PATH)/include -I$(CUDA_PATH)/include -L$(VLI_PATH)/lib -L$(CUDA_PATH)/lib -lvli_gpu -lcuda -lcudart -o bin/high_t__vli_cpu-version
vlicpu:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) -I$(BOOST_PATH)/include minimal_high_t_expansion.cpp -DUSE_VLI_INTEGERS_CPU -I$(VLI_PATH)/include -I$(CUDA_PATH)/include -L$(VLI_PATH)/lib -L$(CUDA_PATH)/lib64 -lvli_gpu -lcuda -lcudart -o bin/high_t__vli_gpu-version
int64:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) minimal_high_t_expansion.cpp -DUSE_INT64_INTEGERS -o bin/high_t__int64-version
gmp:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) minimal_high_t_expansion.cpp -DUSE_GMP_INTEGERS -I$(GMP_PATH)/include -L$(GMP_PATH)/lib -lgmpxx -lgmp -o bin/high_t__gmp-version

run:
	.bin/high_t__vli-version 4 4
	.bin/high_t__gmp-version 4 4

.PHONY: clean
clean:
	-rm bin/*
