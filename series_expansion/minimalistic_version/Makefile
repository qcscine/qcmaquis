#all: int64 gmp vlicpu vligpu
#CC=g++
CC=llvm-g++
CXXFLAGS=-Wall -O3 -mtune=native -funroll-loops -fopenmp
COMMON_DEFS= -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DVLI_MAX_ORDER=4 -DVECTOR_MAX_ORDER=512 


all:  vlicpu gmp

vligpu:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) -I/usr/local/cuda/include -I/opt/boost/include minimal_high_t_expansion_using_dense_vectors.cpp -DUSE_VLI_INTEGERS_GPU -I../../large_integer/VLI/ -L/../../large_integer/VLI/ -L/usr/local/cuda/lib  -o -lcuda -lcudart -lkernels_gpu  bin/high_t__vli_gpu-version
vlicpu:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) -I/opt/boost/include minimal_high_t_expansion_using_dense_vectors.cpp -DUSE_VLI_INTEGERS_CPU -I../../large_integer/VLI/ -o bin/high_t__vli_cpu-version
int64:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) minimal_high_t_expansion_using_dense_vectors.cpp -DUSE_INT64_INTEGERS -o bin/high_t__int64-version
gmp:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) minimal_high_t_expansion_using_dense_vectors.cpp -DUSE_GMP_INTEGERS -lgmp -o bin/high_t__gmp-version

run:
	.bin/high_t__vli-version 4 4
	.bin/high_t__gmp-version 4 4

.PHONY: clean
clean:
	-rm bin/*
