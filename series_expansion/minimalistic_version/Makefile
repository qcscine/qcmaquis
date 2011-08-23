#all: int64-version gmp-version vli-version
all:  vlicpu gmp

vligpu:
	llvm-g++ -I/opt/boost/include  -O3  -g -tune=native -funroll-loops minimal_high_t_expansion_using_dense_vectors_tim.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_VLI_INTEGERS_GPU -I../../large_integer/VLI/ -o bin/high_t__vli-version
	mv bin/high_t__vli-version vlicpu_ex
vlicpu:
	llvm-g++ -I/opt/boost/include  -O3  -g -tune=native -funroll-loops minimal_high_t_expansion_using_dense_vectors_tim.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_VLI_INTEGERS_CPU -I../../large_integer/VLI/ -o bin/high_t__vli-version
	mv bin/high_t__vli-version vlicpu_ex
gmp:
	g++ -Wall -O3 -tune=native -funroll-loops minimal_high_t_expansion_using_dense_vectors_tim.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_GMP_INTEGERS -lgmp -o bin/high_t__gmp-version
	mv bin/high_t__gmp-version gmp_ex

run:
	./gmp_ex 4 4
	./vli_ex 4 4
.PHONY: clean
clean:
	-rm bin/*
