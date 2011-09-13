#all: int64-version gmp-version vli-version
CC=g++
#CC=llvm-g++

BOOST_PATH=/opt/boost
VLI_PATH=/home/hehn/vli_install
CUDA_PATH=/usr/local/cuda

CXXFLAGS=-Wall -O3 -march=native -funroll-loops
#COMMON_DEFS= -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG 
COMMON_DEFS= -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DVLI_MAX_ORDER=4 -DVECTOR_MAX_ORDER=512 

all: vlicpu gmp

vligpu:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) -I$(BOOST_PATH)/include minimal_high_t_expansion_using_dense_vectors.cpp -DUSE_VLI_INTEGERS_GPU -I$(VLI_PATH)/include -I$(CUDA_PATH)/include -L$(VLI_PATH)/lib -lvli_gpu -o bin/high_t__vli-version
vlicpu:
	$(CC) $(CXXFLAGS) $(COMMON_DEFS) -I$(BOOST_PATH)/include minimal_high_t_expansion_using_dense_vectors.cpp -DUSE_VLI_INTEGERS_CPU -I$(VLI_PATH)/include -L$(VLI_PATH)/lib -lvli_gpu -o bin/high_t__vli-version
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
