#include "utils/testing.hpp"

TEST_CASE( "Matrix QR factorization is computed", "[qr]" )
{
    matrix<double> A(TEST_M,TEST_M);
    matrix<double> Q(TEST_M,TEST_M);
    matrix<double> R(TEST_M,TEST_M);
    matrix<double> C(TEST_M,TEST_M);

    generate(A);
    qr(A, Q, R);
    gemm(Q, R, C);

    REQUIRE((C == A));
}
