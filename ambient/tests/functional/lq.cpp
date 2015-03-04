#include "utils/testing.hpp"

TEST_CASE( "Matrix LQ factorization is computed", "[lq]" )
{
    matrix<double> L(TEST_M,TEST_N);
    matrix<double> Q(TEST_M,TEST_N);
    matrix<double> A(TEST_M,TEST_N);
    matrix<double> C(TEST_M,TEST_N);

    generate(A);
    lq(A, L, Q);
    gemm(L, Q, C);

    REQUIRE((C == A));
}
