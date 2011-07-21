all: int64-version gmp-version vli-version

#int64-version:
#	g++ -Wall -O3 -march=native -funroll-loops minimal_high_t_expansion_using_sparse_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version
#
#gmp-version:
#	g++ -Wall -O3 -march=native -funroll-loops minimal_high_t_expansion_using_sparse_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_GMP_INTEGERS -lgmp -o bin/high_t__gmp-version

vli-version:
	g++ -Wall -O5 -march=native -funroll-loops minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_VLI_INTEGERS -I../../large_integer/VLI/ -o bin/high_t__vli-version
#	g++ -Wall -O0 -g minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version

int64-version:
	g++ -Wall -O5 -march=native -funroll-loops minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version
#	g++ -Wall -O0 -g minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version

gmp-version:
	g++ -Wall -O5 -march=native -funroll-loops minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_GMP_INTEGERS -lgmp -o bin/high_t__gmp-version

.PHONY: clean
clean:
	-rm bin/*
