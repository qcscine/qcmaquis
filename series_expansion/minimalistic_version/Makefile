all: int64-version gmp-version

#int64-version:
#	g++ -Wall -O2 minimal_high_t_expansion_using_sparse_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version
#
#gmp-version:
#	g++ -Wall -O2 minimal_high_t_expansion_using_sparse_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_GMP_INTEGERS -lgmp -o bin/high_t__gmp-version

int64-version:
	g++ -Wall -O2 minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version
#	g++ -Wall -O0 -g minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_INT64_INTEGERS -o bin/high_t__int64-version

gmp-version:
	g++ -Wall -O2 minimal_high_t_expansion_using_dense_vectors.cpp -DPOLYNOMIAL_MAX_ORDER=21 -DNDEBUG -DUSE_GMP_INTEGERS -lgmp -o bin/high_t__gmp-version

.PHONY: clean
clean:
	-rm bin/*
